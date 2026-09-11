from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
IN_FILE = PROJECT_DIR / "outputs/reports/unknown_chembl_screen/top100_putative_new_hits.csv"
OUT_DIR = PROJECT_DIR / "outputs/reports/unknown_chembl_screen"
TOP20_DIVERSE = OUT_DIR / "top20_diverse_putative_new_hits.csv"
SUMMARY_FILE = OUT_DIR / "top_hit_prioritization_summary.csv"


########## 1. load ##########
df = pd.read_csv(IN_FILE)
df["admet_lite_pass"] = (
    (df["QED"] >= 0.55)
    & (df["MW"].between(200, 450))
    & (df["LogP"].between(1.0, 4.5))
    & (df["TPSA"].between(30, 100))
    & (df["RotBonds"] <= 8)
    & (df["PAINS_count"] == 0)
    & (df["Lipinski_violations"] == 0)
    & (df["Veber_pass"])
)
df["priority_score"] = (
    0.45 * df["trpv1_binding_score"]
    + 0.30 * df["v1_selectivity_score"]
    + 0.20 * df["QED"]
    + 0.05 * df["novel_scaffold_to_trp_train"].astype(int)
)
df["tier"] = "B"
df.loc[
    (df["trpv1_binding_score"] >= 0.85)
    & (df["v1_selectivity_score"] >= 0.85)
    & df["admet_lite_pass"],
    "tier",
] = "A"
df.loc[
    (df["trpv1_binding_score"] < 0.75)
    | (df["v1_selectivity_score"] < 0.75)
    | (~df["admet_lite_pass"]),
    "tier",
] = "C"


########## 2. diverse selection ##########
ranked = df.sort_values(
    ["tier", "priority_score", "trpv1_binding_score", "v1_selectivity_score"],
    ascending=[True, False, False, False],
).reset_index(drop=True)
diverse_rows = []
seen_scaffolds = set()
for _, row in ranked.iterrows():
    scaffold = row["scaffold"]
    if scaffold in seen_scaffolds:
        continue
    diverse_rows.append(row)
    seen_scaffolds.add(scaffold)
    if len(diverse_rows) >= 20:
        break
top20 = pd.DataFrame(diverse_rows)
top20.to_csv(TOP20_DIVERSE, index=False)


########## 3. plots ##########
colors = {"A": "#2a9d8f", "B": "#457b9d", "C": "#c1666b"}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for tier, sub in df.groupby("tier"):
    axes[0].scatter(
        sub["trpv1_binding_score"],
        sub["v1_selectivity_score"],
        s=35 + 60 * sub["QED"],
        alpha=0.75,
        label=f"Tier {tier}",
        color=colors[tier],
    )
axes[0].set_xlabel("TRPV1 binding score")
axes[0].set_ylabel("V1 selectivity score")
axes[0].legend(frameon=False)
axes[0].set_title("Model score prioritization")
axes[1].scatter(df["LogP"], df["QED"], c=df["TPSA"], cmap="viridis", s=55, alpha=0.8)
axes[1].set_xlabel("LogP")
axes[1].set_ylabel("QED")
axes[1].set_title("Drug-likeness space")
fig.tight_layout()
fig.savefig(OUT_DIR / "unknown_hit_prioritization.png", dpi=300)
plt.close(fig)

mols = [Chem.MolFromSmiles(s) for s in top20["canonical_smiles"]]
legends = [
    f"{chembl_id}\\nT{tier} P={priority:.2f}"
    for chembl_id, tier, priority in zip(
        top20["molecule_chembl_id"],
        top20["tier"],
        top20["priority_score"],
    )
]
img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(320, 220), legends=legends)
img.save(str(OUT_DIR / "top20_diverse_putative_new_hit_structures.png"))


########## 4. save ##########
summary = (
    df.groupby(["tier", "admet_lite_pass"])
    .size()
    .reset_index(name="n")
    .sort_values(["tier", "admet_lite_pass"])
)
summary.to_csv(SUMMARY_FILE, index=False)
print(f"Saved: {TOP20_DIVERSE}")
print(f"Saved: {SUMMARY_FILE}")

