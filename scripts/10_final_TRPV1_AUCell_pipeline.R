suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(hdf5r)
  library(Matrix)
  library(ggplot2)
})

# Final TRPV1 AUCell pipeline
# Scope:
#   1. Use only target diseases ASD, FTD, HD, MS, LBD, and PD.
#   2. Analyze only TRPV1-downregulated disease-subgroup comparisons:
#      p-value < 0.05 and logFC < -0.2.
#   3. Score each cell from the selected disease subgroup and its matching CTRL
#      subgroup for protein aggregation/proteostasis, mitochondrial stress, chronic
#      inflammation, and related pathway gene sets.
#   4. Save tables and figures under project outputs.

project_dir <- "/home/nizhu/Projects/06.TRPV1_vs_screening"
root <- file.path(project_dir, "outputs", "geo_trpv1")
h5ad_path <- file.path(project_dir, "human_cortex_restricted.h5ad")
de_path <- file.path(root, "TRPV1_DE_vs_CT.csv")
pb_path <- file.path(root, "TRPV1_donor_pseudobulk_summary.csv")

table_dir <- root
fig_dir <- root
logs_dir <- file.path(project_dir, "logs")
bpcells_dir <- file.path(root, "bpcells_h5ad_colmajor")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(de_path) && file.exists(file.path(root, "TRPV1_DE_vs_CT.csv"))) {
  file.copy(file.path(root, "TRPV1_DE_vs_CT.csv"), de_path)
}
if (!file.exists(pb_path) && file.exists(file.path(root, "TRPV1_donor_pseudobulk_summary.csv"))) {
  file.copy(file.path(root, "TRPV1_donor_pseudobulk_summary.csv"), pb_path)
}

target_diseases <- c("ASD", "FTD", "HD", "MS", "LBD", "PD")
chunk_size <- as.integer(Sys.getenv("AUCELL_CHUNK_SIZE", "5000"))
auc_max_rank <- as.integer(Sys.getenv("AUCELL_MAX_RANK", "1000"))

message("Root: ", root)
message("Chunk size: ", chunk_size)
message("AUC max rank: ", auc_max_rank)

read_h5ad_obs_column <- function(h5, column) {
  node <- h5[["obs"]][[column]]
  if (inherits(node, "H5Group")) {
    if ("codes" %in% names(node) && "categories" %in% names(node)) {
      codes <- node[["codes"]][]
      cats <- node[["categories"]][]
      out <- rep(NA_character_, length(codes))
      keep <- !is.na(codes) & codes >= 0
      out[keep] <- as.character(cats[codes[keep] + 1L])
      return(out)
    }
    if ("values" %in% names(node)) {
      return(as.character(node[["values"]][]))
    }
  }
  as.character(node[])
}

sanitize_id <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

score_file_for <- function(comparison_id) {
  file.path(table_dir, paste0("AUCell_scores_", sanitize_id(comparison_id), ".csv"))
}

theme_of_gene_set <- function(gene_set) {
  ifelse(grepl("protein|unfolded", gene_set), "Protein aggregation / proteostasis",
    ifelse(grepl("mitochondrial|oxidative", gene_set), "Mitochondrial stress",
      "Chronic inflammation"))
}

save_png <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300)
}

gene_sets <- list(
  protein_aggregation_proteostasis = c(
    "HSPA1A","HSPA1B","HSPA5","HSPA8","HSP90AA1","HSP90AB1","HSPB1","DNAJA1",
    "DNAJB1","DNAJB6","HSPH1","HSF1","BAG3","STUB1","VCP","SQSTM1","UBC",
    "UBB","PSMA1","PSMA2","PSMA3","PSMB1","PSMB5","PSMD1","PSMD2","ATG5",
    "ATG7","BECN1","MAP1LC3B","LAMP2"
  ),
  unfolded_protein_response_er_stress = c(
    "HSPA5","HSP90B1","CALR","CANX","PDIA3","PDIA4","XBP1","ERN1","EIF2AK3",
    "ATF4","ATF6","DDIT3","DNAJC3","HERPUD1","EDEM1","SEL1L","DERL1","SYVN1",
    "PPP1R15A","TRIB3"
  ),
  mitochondrial_stress_uprmt = c(
    "HSPD1","HSPE1","HSPA9","LONP1","CLPP","CLPX","YME1L1","OMA1","ATF5",
    "DDIT3","CHCHD10","SOD2","PRDX3","TXN2","GPX4","PINK1","PRKN","BNIP3",
    "BNIP3L","FUNDC1","TOMM20","TIMM23"
  ),
  oxidative_phosphorylation_ros = c(
    "NDUFA1","NDUFA2","NDUFA4","NDUFB3","NDUFB8","NDUFS1","NDUFS2","SDHA",
    "SDHB","UQCRC1","UQCRC2","COX4I1","COX5A","COX6A1","ATP5F1A","ATP5F1B",
    "ATP5MC1","SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","PRDX3","TXN","TXN2",
    "NFE2L2","HMOX1","NQO1"
  ),
  chronic_inflammation_nfkb_cytokine = c(
    "NFKB1","NFKB2","RELA","RELB","CHUK","IKBKB","NFKBIA","TNF","TNFRSF1A",
    "IL1B","IL1A","IL6","IL18","CXCL8","CXCL10","CCL2","CCL3","CCL4","CCL5",
    "PTGS2","ICAM1","VCAM1","STAT1","STAT3","SOCS3","JUN","FOS"
  ),
  inflammasome_complement_microglia = c(
    "NLRP3","PYCARD","CASP1","IL1B","IL18","GSDMD","AIM2","NLRC4","C1QA",
    "C1QB","C1QC","C3","C4A","C4B","CFB","ITGAM","TYROBP","TREM2","AIF1",
    "LST1","FCGR3A","CD68","CTSS","HLA-DRA","HLA-DRB1"
  ),
  interferon_antiviral_inflammation = c(
    "IFIT1","IFIT2","IFIT3","ISG15","MX1","MX2","OAS1","OAS2","OAS3","IFI6",
    "IFI27","IFI44","IFI44L","IRF1","IRF7","IRF9","STAT1","STAT2","CXCL10",
    "GBP1","GBP2","RSAD2","DDX58","IFIH1"
  )
)

de <- read.csv(de_path, check.names = FALSE)
pb <- read.csv(pb_path, check.names = FALSE)
de$disease <- sub("-CT$", "", de[["contrast"]])

hits <- de[
  de[["Gene"]] == "TRPV1" &
    de$disease %in% target_diseases &
    as.numeric(de[["p-value"]]) < 0.05 &
    as.numeric(de[["logFC"]]) < -0.2,
]
hits$comparison_id <- paste(hits$disease, hits[["Cell type"]], sep = "__")
hits$score_file <- score_file_for(hits$comparison_id)
write.csv(hits, file.path(table_dir, "TRPV1_downregulated_target_subgroups_ASD_FTD_HD_MS_LBD_PD.csv"), row.names = FALSE)

disease_availability <- data.frame(
  disease = target_diseases,
  present_in_DE_contrast = target_diseases %in% unique(de$disease),
  has_downregulated_hit = target_diseases %in% unique(hits$disease)
)
write.csv(disease_availability, file.path(table_dir, "target_disease_availability.csv"), row.names = FALSE)

gene_set_df <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  data.frame(gene_set = nm, gene = gene_sets[[nm]], stringsAsFactors = FALSE)
}))
write.csv(gene_set_df, file.path(table_dir, "AUCell_gene_sets_used.csv"), row.names = FALSE)

if (!dir.exists(bpcells_dir)) {
  message("Creating BPCells column-major matrix directory: ", bpcells_dir)
  Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
  mat0 <- open_matrix_anndata_hdf5(h5ad_path, group = "X")
  transpose_storage_order(
    mat0,
    outdir = bpcells_dir,
    tmpdir = file.path(root, "tmp_transpose"),
    sort_bytes = 2^30
  )
}

mat <- open_matrix_dir(bpcells_dir)
genes_present <- rownames(mat)
gene_sets <- lapply(gene_sets, function(gs) intersect(unique(gs), genes_present))
gene_set_sizes <- data.frame(gene_set = names(gene_sets), n_genes_present = lengths(gene_sets))
write.csv(gene_set_sizes, file.path(table_dir, "AUCell_gene_set_sizes_present.csv"), row.names = FALSE)
gene_sets <- gene_sets[lengths(gene_sets) >= 5]

score_method <- data.frame(
  method = "chunked_AUCell_formula",
  description = paste0(
    "For each cell, non-zero genes are ranked by expression descending. ",
    "For each gene set, the score is the normalized recovery AUC of member genes ",
    "within the top aucMaxRank genes: sum(aucMaxRank-rank+1)/(aucMaxRank*min(gene_set_size, aucMaxRank))."
  ),
  aucMaxRank = auc_max_rank,
  chunk_size = chunk_size
)
write.csv(score_method, file.path(table_dir, "AUCell_scoring_method.csv"), row.names = FALSE)

h5 <- H5File$new(h5ad_path, mode = "r")
obs <- data.frame(
  cell = read_h5ad_obs_column(h5, "_index"),
  condition = read_h5ad_obs_column(h5, "condition"),
  sub_cell_type = read_h5ad_obs_column(h5, "sub_cell_type"),
  donor = read_h5ad_obs_column(h5, "donor"),
  stringsAsFactors = FALSE
)
h5$close_all()

write.csv(
  aggregate(n_cells ~ sub_cell_type + condition, pb, sum),
  file.path(table_dir, "input_pseudobulk_cell_counts_by_subtype_condition.csv"),
  row.names = FALSE
)

expected <- do.call(rbind, lapply(seq_len(nrow(hits)), function(i) {
  disease <- hits$disease[i]
  subtype <- hits[["Cell type"]][i]
  expected_cells <- sum(pb$n_cells[pb$sub_cell_type == subtype & pb$condition %in% c("CT", disease)])
  file_rows <- if (file.exists(hits$score_file[i])) max(0L, length(readLines(hits$score_file[i])) - 1L) else NA_integer_
  data.frame(
    disease = disease,
    sub_cell_type = subtype,
    comparison_id = hits$comparison_id[i],
    file = hits$score_file[i],
    expected_cells = expected_cells,
    file_rows = file_rows,
    complete = !is.na(file_rows) && as.integer(expected_cells) == as.integer(file_rows),
    stringsAsFactors = FALSE
  )
}))
write.csv(expected, file.path(table_dir, "resume_audit_target_expected_vs_existing.csv"), row.names = FALSE)

score_one_chunk <- function(cell_idx, hit_row) {
  expr <- as(mat[, cell_idx, drop = FALSE], "dgCMatrix")
  gene_index_to_row <- seq_len(nrow(expr))
  names(gene_index_to_row) <- rownames(expr)
  local_gene_set_indices <- lapply(gene_sets, function(gs) unname(gene_index_to_row[gs]))
  local_gene_set_indices <- lapply(local_gene_set_indices, function(x) x[!is.na(x)])
  score_mat <- matrix(
    0,
    nrow = length(local_gene_set_indices),
    ncol = ncol(expr),
    dimnames = list(names(local_gene_set_indices), colnames(expr))
  )

  for (col_i in seq_len(ncol(expr))) {
    start <- expr@p[col_i] + 1L
    end <- expr@p[col_i + 1L]
    if (start > end) next
    rows <- expr@i[start:end] + 1L
    vals <- expr@x[start:end]
    ord <- order(vals, decreasing = TRUE)
    rows <- rows[ord]
    top_n <- min(length(rows), auc_max_rank)
    if (top_n == 0L) next
    top_rows <- rows[seq_len(top_n)]
    ranks <- seq_len(top_n)
    names(ranks) <- as.character(top_rows)

    for (gs_i in seq_along(local_gene_set_indices)) {
      gs_rows <- local_gene_set_indices[[gs_i]]
      hit_ranks <- ranks[as.character(gs_rows)]
      hit_ranks <- hit_ranks[!is.na(hit_ranks)]
      if (length(hit_ranks) == 0L) next
      denom <- auc_max_rank * min(length(gs_rows), auc_max_rank)
      score_mat[gs_i, col_i] <- sum(auc_max_rank - hit_ranks + 1) / denom
    }
  }

  scores <- as.data.frame(t(score_mat), check.names = FALSE)
  scores$cell <- colnames(expr)
  scores$sub_cell_type <- hit_row[["Cell type"]]
  scores$disease <- hit_row[["disease"]]
  scores$contrast <- hit_row[["contrast"]]
  scores$logFC_TRPV1 <- as.numeric(hit_row[["logFC"]])
  scores$p_value_TRPV1 <- as.numeric(hit_row[["p-value"]])
  scores$p_adj_TRPV1 <- as.numeric(hit_row[["p-value adj."]])
  scores$condition <- obs$condition[cell_idx]
  scores$donor <- obs$donor[cell_idx]
  scores
}

for (i in seq_len(nrow(hits))) {
  if (isTRUE(expected$complete[i])) {
    message("Reusing complete score file: ", basename(hits$score_file[i]))
    next
  }

  hit <- hits[i, , drop = FALSE]
  subtype <- hit[["Cell type"]]
  disease <- hit[["disease"]]
  comparison_id <- hit[["comparison_id"]]
  score_path <- hit[["score_file"]]
  message("Scoring ", comparison_id)

  keep <- obs$sub_cell_type == subtype & obs$condition %in% c("CT", disease)
  idx <- which(keep)
  chunks <- split(idx, ceiling(seq_along(idx) / chunk_size))
  if (file.exists(score_path)) file.remove(score_path)

  for (j in seq_along(chunks)) {
    message("  chunk ", j, "/", length(chunks), " cells=", length(chunks[[j]]))
    chunk_scores <- score_one_chunk(chunks[[j]], hit)
    write.table(
      chunk_scores,
      file = score_path,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(score_path),
      append = file.exists(score_path),
      quote = TRUE
    )
    rm(chunk_scores)
    gc()
  }
}

score_files <- hits$score_file
writeLines(score_files, file.path(table_dir, "AUCell_target_score_files.txt"))

summary_rows <- list()
donor_rows <- list()
for (i in seq_len(nrow(hits))) {
  scores <- read.csv(hits$score_file[i], check.names = FALSE)
  meta_cols <- c("cell", "sub_cell_type", "disease", "contrast", "logFC_TRPV1",
                 "p_value_TRPV1", "p_adj_TRPV1", "condition", "donor")
  score_gene_sets <- setdiff(colnames(scores), meta_cols)
  disease <- hits$disease[i]
  subtype <- hits[["Cell type"]][i]
  comparison_id <- hits$comparison_id[i]

  for (gs in score_gene_sets) {
    disease_scores <- scores[scores$condition == disease, gs]
    ctrl_scores <- scores[scores$condition == "CT", gs]
    wt <- suppressWarnings(tryCatch(wilcox.test(disease_scores, ctrl_scores), error = function(e) NULL))
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      comparison_id = comparison_id,
      disease = disease,
      sub_cell_type = subtype,
      contrast = hits$contrast[i],
      TRPV1_logFC = as.numeric(hits$logFC[i]),
      TRPV1_p_value = as.numeric(hits[["p-value"]][i]),
      TRPV1_p_adj = as.numeric(hits[["p-value adj."]][i]),
      gene_set = gs,
      n_disease_cells = sum(scores$condition == disease),
      n_ctrl_cells = sum(scores$condition == "CT"),
      mean_disease = mean(disease_scores, na.rm = TRUE),
      mean_ctrl = mean(ctrl_scores, na.rm = TRUE),
      median_disease = median(disease_scores, na.rm = TRUE),
      median_ctrl = median(ctrl_scores, na.rm = TRUE),
      delta_mean_disease_minus_ctrl = mean(disease_scores, na.rm = TRUE) - mean(ctrl_scores, na.rm = TRUE),
      wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value,
      stringsAsFactors = FALSE
    )
  }

  long_for_donor <- do.call(rbind, lapply(score_gene_sets, function(gs) {
    data.frame(
      comparison_id = comparison_id,
      disease = disease,
      sub_cell_type = subtype,
      condition = scores$condition,
      donor = scores$donor,
      gene_set = gs,
      score = scores[[gs]],
      stringsAsFactors = FALSE
    )
  }))
  donor_rows[[length(donor_rows) + 1L]] <- aggregate(
    score ~ comparison_id + disease + sub_cell_type + condition + donor + gene_set,
    long_for_donor,
    mean
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_df$wilcox_p_adj_BH <- p.adjust(summary_df$wilcox_p, method = "BH")
summary_df$disease_higher <- summary_df$delta_mean_disease_minus_ctrl > 0
summary_df$theme <- theme_of_gene_set(summary_df$gene_set)
summary_df$significant_disease_higher <- summary_df$disease_higher & summary_df$wilcox_p_adj_BH < 0.05
write.csv(summary_df, file.path(table_dir, "AUCell_target_disease_vs_CTRL_summary.csv"), row.names = FALSE)
write.csv(summary_df[summary_df$disease_higher, ], file.path(table_dir, "AUCell_target_disease_higher_scores.csv"), row.names = FALSE)
write.csv(
  summary_df[summary_df$significant_disease_higher, ],
  file.path(table_dir, "AUCell_target_disease_higher_scores_BH0.05.csv"),
  row.names = FALSE
)

comparison_summary <- do.call(rbind, lapply(split(summary_df, summary_df$comparison_id), function(x) {
  x_ord <- x[order(x$wilcox_p_adj_BH, -x$delta_mean_disease_minus_ctrl), ]
  sig <- x_ord[x_ord$significant_disease_higher, ]
  data.frame(
    comparison_id = x$comparison_id[1],
    disease = x$disease[1],
    sub_cell_type = x$sub_cell_type[1],
    TRPV1_logFC = x$TRPV1_logFC[1],
    TRPV1_p_value = x$TRPV1_p_value[1],
    n_gene_sets_tested = nrow(x),
    n_disease_higher = sum(x$disease_higher, na.rm = TRUE),
    n_significant_disease_higher_BH0.05 = nrow(sig),
    top_disease_higher_gene_set = if (nrow(sig) == 0) NA_character_ else sig$gene_set[1],
    top_delta_mean = if (nrow(sig) == 0) NA_real_ else sig$delta_mean_disease_minus_ctrl[1],
    top_BH_p = if (nrow(sig) == 0) NA_real_ else sig$wilcox_p_adj_BH[1],
    significant_gene_sets = if (nrow(sig) == 0) "" else paste(sig$gene_set, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
comparison_summary <- comparison_summary[order(comparison_summary$disease, comparison_summary$sub_cell_type), ]
write.csv(comparison_summary, file.path(table_dir, "AUCell_final_by_TRPV1_downregulated_subgroup.csv"), row.names = FALSE)

theme_summary <- aggregate(
  cbind(delta_mean_disease_minus_ctrl, disease_higher, significant_disease_higher) ~ disease + sub_cell_type + theme,
  data = summary_df,
  FUN = function(z) if (is.logical(z)) sum(z, na.rm = TRUE) else mean(z, na.rm = TRUE)
)
names(theme_summary)[names(theme_summary) == "delta_mean_disease_minus_ctrl"] <- "mean_delta_across_gene_sets"
names(theme_summary)[names(theme_summary) == "disease_higher"] <- "n_gene_sets_disease_higher"
names(theme_summary)[names(theme_summary) == "significant_disease_higher"] <- "n_gene_sets_significant_disease_higher_BH0.05"
write.csv(theme_summary, file.path(table_dir, "AUCell_final_theme_summary.csv"), row.names = FALSE)

donor_df <- do.call(rbind, donor_rows)
write.csv(donor_df, file.path(table_dir, "AUCell_final_donor_mean_scores.csv"), row.names = FALSE)

donor_stats <- do.call(rbind, lapply(split(donor_df, paste(donor_df$comparison_id, donor_df$gene_set, sep = "||")), function(x) {
  disease <- x$disease[1]
  d <- x$score[x$condition == disease]
  c <- x$score[x$condition == "CT"]
  wt <- suppressWarnings(tryCatch(wilcox.test(d, c), error = function(e) NULL))
  data.frame(
    comparison_id = x$comparison_id[1],
    disease = disease,
    sub_cell_type = x$sub_cell_type[1],
    gene_set = x$gene_set[1],
    n_disease_donors = length(unique(x$donor[x$condition == disease])),
    n_ctrl_donors = length(unique(x$donor[x$condition == "CT"])),
    mean_disease_donor_score = mean(d, na.rm = TRUE),
    mean_ctrl_donor_score = mean(c, na.rm = TRUE),
    delta_mean_disease_minus_ctrl = mean(d, na.rm = TRUE) - mean(c, na.rm = TRUE),
    wilcox_p_donor = if (is.null(wt)) NA_real_ else wt$p.value,
    stringsAsFactors = FALSE
  )
}))
donor_stats$wilcox_p_donor_BH <- p.adjust(donor_stats$wilcox_p_donor, method = "BH")
donor_stats$disease_higher_donor_mean <- donor_stats$delta_mean_disease_minus_ctrl > 0
write.csv(donor_stats, file.path(table_dir, "AUCell_final_donor_level_disease_vs_CTRL.csv"), row.names = FALSE)

summary_df$neglog10_BH <- -log10(pmax(summary_df$wilcox_p_adj_BH, .Machine$double.xmin))
summary_df$comparison_label <- paste(summary_df$disease, summary_df$sub_cell_type, sep = " | ")

p0 <- ggplot(summary_df, aes(
  x = gene_set,
  y = reorder(comparison_label, disease),
  fill = delta_mean_disease_minus_ctrl
)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030", midpoint = 0) +
  labs(
    title = "AUCell pathway score difference in TRPV1-downregulated subgroups",
    subtitle = "Rows are disease-subgroup comparisons; columns are target pathway gene sets; red indicates disease group higher than CTRL",
    x = "Target pathway gene set",
    y = "Disease | TRPV1-downregulated cell subgroup",
    fill = "Mean score\nDisease - CTRL"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_png(p0, "AUCell_target_delta_heatmap", 11, 8)

p1 <- ggplot(summary_df, aes(x = gene_set, y = reorder(comparison_label, disease),
                             size = neglog10_BH, color = delta_mean_disease_minus_ctrl)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(low = "#2f6f9f", mid = "grey92", high = "#b43c3c", midpoint = 0) +
  scale_size_continuous(range = c(1.5, 7)) +
  labs(
    title = "AUCell pathway score dot plot for TRPV1-downregulated subgroups",
    subtitle = "Color shows mean disease-vs-CTRL score difference; point size shows BH-adjusted significance",
    x = "Target pathway gene set",
    y = "Disease | TRPV1-downregulated cell subgroup",
    color = "Mean score\nDisease - CTRL",
    size = "-log10(BH)"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
save_png(p1, "AUCell_tutorial_style_dotplot_all_target_comparisons", 11, 8)

p2 <- ggplot(subset(summary_df, significant_disease_higher), aes(
  x = gene_set,
  y = reorder(comparison_label, delta_mean_disease_minus_ctrl),
  fill = delta_mean_disease_minus_ctrl
)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#f6d4bf", high = "#9f2f2f") +
  labs(
    title = "Significant disease-higher AUCell pathway scores",
    subtitle = "Only comparisons with disease mean > CTRL mean and BH-adjusted P < 0.05 are displayed",
    x = "Target pathway gene set",
    y = "Disease | TRPV1-downregulated cell subgroup",
    fill = "Mean score\nDisease - CTRL"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_png(p2, "AUCell_tutorial_style_significant_up_heatmap", 10, 7.5)

theme_summary$comparison_label <- paste(theme_summary$disease, theme_summary$sub_cell_type, sep = " | ")
p3 <- ggplot(theme_summary, aes(x = theme, y = comparison_label, fill = mean_delta_across_gene_sets)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#c53030", midpoint = 0) +
  labs(
    title = "Mean AUCell score difference by pathway theme",
    subtitle = "Each tile averages related gene sets within a pathway theme; red indicates disease group higher than CTRL",
    x = "Pathway theme",
    y = "Disease | TRPV1-downregulated cell subgroup",
    fill = "Mean delta\nDisease - CTRL"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
save_png(p3, "AUCell_final_theme_delta_heatmap", 9, 7)

comparison_summary$label <- factor(paste(comparison_summary$disease, comparison_summary$sub_cell_type, sep = " | "),
                                   levels = rev(paste(comparison_summary$disease, comparison_summary$sub_cell_type, sep = " | ")))
p4 <- ggplot(comparison_summary, aes(x = n_significant_disease_higher_BH0.05, y = label, fill = disease)) +
  geom_col() +
  labs(
    title = "Number of target pathway gene sets with higher disease-group AUCell scores",
    subtitle = "Only TRPV1-downregulated disease-subgroup comparisons are shown; counts use BH-adjusted P < 0.05 and disease mean > CTRL mean",
    x = "Number of significant disease-higher gene sets",
    y = "Disease | TRPV1-downregulated cell subgroup",
    fill = "Disease"
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank())
save_png(p4, "AUCell_final_significant_gene_set_counts", 9, 7)

top40 <- summary_df[order(summary_df$wilcox_p_adj_BH, -abs(summary_df$delta_mean_disease_minus_ctrl)), ]
top40 <- head(top40, 40)
top40$label <- factor(paste(top40$comparison_label, top40$gene_set, sep = " | "),
                      levels = rev(paste(top40$comparison_label, top40$gene_set, sep = " | ")))
p5 <- ggplot(top40, aes(x = delta_mean_disease_minus_ctrl, y = label, fill = disease_higher)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.25) +
  scale_fill_manual(values = c("FALSE" = "#2b6cb0", "TRUE" = "#c53030"), labels = c("CTRL higher", "Disease higher")) +
  labs(
    title = "Top AUCell pathway score differences",
    subtitle = "Top 40 disease-subgroup-gene-set comparisons ranked by adjusted significance and effect size",
    x = "Mean score difference: disease group - CTRL",
    y = "Disease | subgroup | gene set",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank())
save_png(p5, "AUCell_target_top40_delta_barplot", 12, 8)

top_rows <- summary_df[summary_df$significant_disease_higher, ]
top_rows <- top_rows[order(top_rows$theme, top_rows$disease, -top_rows$delta_mean_disease_minus_ctrl), ]
top_rows <- do.call(rbind, lapply(split(top_rows, top_rows$theme), head, 6))
write.csv(top_rows, file.path(table_dir, "AUCell_plotting_top_significant_for_violin.csv"), row.names = FALSE)

dist_df <- do.call(rbind, lapply(seq_len(nrow(top_rows)), function(i) {
  row <- top_rows[i, ]
  x <- read.csv(score_file_for(row$comparison_id), check.names = FALSE)
  data.frame(
    comparison_id = row$comparison_id,
    disease = row$disease,
    sub_cell_type = row$sub_cell_type,
    gene_set = row$gene_set,
    theme = row$theme,
    condition = x$condition,
    score = x[[row$gene_set]],
    stringsAsFactors = FALSE
  )
}))
dist_df$panel <- paste(dist_df$comparison_id, dist_df$gene_set, sep = "\n")
write.csv(dist_df, file.path(table_dir, "AUCell_plotting_violin_scores_top_significant.csv"), row.names = FALSE)

p6 <- ggplot(dist_df, aes(x = condition, y = score, fill = condition)) +
  geom_violin(scale = "width", linewidth = 0.15, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.size = 0.05, linewidth = 0.15, alpha = 0.65) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("CT" = "#6b8fb3", "ASD" = "#c65f4b", "FTD" = "#c65f4b", "HD" = "#c65f4b", "MS" = "#c65f4b")) +
  labs(
    title = "Single-cell AUCell score distributions for top significant disease-higher pathways",
    subtitle = "Each panel compares cells from one TRPV1-downregulated subgroup between disease and CTRL",
    x = "Condition",
    y = "AUCell-style pathway score"
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none", strip.text = element_text(size = 7), panel.grid.minor = element_blank())
save_png(p6, "AUCell_tutorial_style_violin_top_significant", 11.5, 9.5)

p7 <- ggplot(theme_summary, aes(x = mean_delta_across_gene_sets, y = reorder(comparison_label, mean_delta_across_gene_sets), fill = theme)) +
  geom_col(position = position_dodge2(width = 0.75, preserve = "single"), width = 0.65) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.25) +
  labs(
    title = "Average disease-vs-CTRL AUCell score difference by pathway theme",
    subtitle = "Bars summarize protein aggregation/proteostasis, mitochondrial stress, and chronic inflammation themes",
    x = "Mean score difference across related gene sets: disease - CTRL",
    y = "Disease | TRPV1-downregulated cell subgroup",
    fill = "Pathway theme"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.y = element_blank())
save_png(p7, "AUCell_tutorial_style_theme_delta_barplot", 11, 8)

unlink(file.path(fig_dir, "*.pdf"))

readme <- c(
  "# TRPV1 downregulated subgroup AUCell scoring",
  "",
  "Scope: ASD, FTD, HD, MS, LBD, PD. Only disease-subgroup comparisons with TRPV1 p-value < 0.05 and logFC < -0.2 were analyzed.",
  "",
  "LBD and PD were requested but are not present in the DE contrast table, so no TRPV1-downregulated subgroups were available for these diseases.",
  "",
  "Primary table: outputs/geo_trpv1/AUCell_target_disease_vs_CTRL_summary.csv",
  "Disease-higher significant results: outputs/geo_trpv1/AUCell_target_disease_higher_scores_BH0.05.csv",
  "Final subgroup interpretation: outputs/geo_trpv1/AUCell_final_by_TRPV1_downregulated_subgroup.csv",
  "Theme summary: outputs/geo_trpv1/AUCell_final_theme_summary.csv",
  "Donor-level check: outputs/geo_trpv1/AUCell_final_donor_level_disease_vs_CTRL.csv",
  "",
  "Figures and tables are stored under outputs/geo_trpv1. Logs are stored under logs."
)
writeLines(readme, file.path(root, "README_AUCell_final_outputs.md"))

message("Final comparisons: ", nrow(hits))
message("Gene set comparisons: ", nrow(summary_df))
message("Significant disease-higher rows: ", sum(summary_df$significant_disease_higher))
message("PNG figures: ", length(list.files(fig_dir, pattern = "\\.png$")))
message("DONE")
