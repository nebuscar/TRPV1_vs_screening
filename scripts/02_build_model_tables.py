from pathlib import Path
import pandas as pd


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
INPUT_FILE = PROJECT_DIR / "outputs/ligands/chembl_trp_activity_raw.csv"
OUT_DIR = PROJECT_DIR / "outputs/ligands"
ACTIVITY_TABLE = OUT_DIR / "trp_activity_model_table.csv"
TRPV1_BINDING_TABLE = OUT_DIR / "trpv1_binding_binary_table.csv"
V1_SELECTIVITY_TABLE = OUT_DIR / "trpv1_vs_nonv1_selectivity_table.csv"


########## 1. load ##########
df = pd.read_csv(INPUT_FILE)
df = df[df["canonical_smiles"].notna()].copy()
df["binder_label"] = pd.to_numeric(df["binder_label"], errors="coerce")
df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")


########## 2. collapse activity ##########
activity = (
    df.sort_values(["canonical_smiles", "query_target", "binder_label", "pchembl_value"], ascending=[True, True, False, False])
    .groupby(["canonical_smiles", "query_target"], as_index=False)
    .agg(
        target_role=("target_role", "first"),
        max_pchembl=("pchembl_value", "max"),
        min_standard_value=("standard_value", "min"),
        binder_label=("binder_label", "max"),
        n_records=("molecule_chembl_id", "count"),
        molecule_chembl_id=("molecule_chembl_id", "first"),
    )
)
activity["v1_label"] = (activity["query_target"] == "TRPV1").astype(int)
OUT_DIR.mkdir(parents=True, exist_ok=True)
activity.to_csv(ACTIVITY_TABLE, index=False)


########## 3. trpv1 binding table ##########
trpv1_binding = activity[
    (activity["query_target"] == "TRPV1") & activity["binder_label"].isin([0, 1])
].copy()
trpv1_binding["label"] = trpv1_binding["binder_label"].astype(int)
trpv1_binding.to_csv(TRPV1_BINDING_TABLE, index=False)


########## 4. v1 selectivity table ##########
active = activity[activity["binder_label"] == 1].copy()
compound_targets = (
    active.groupby("canonical_smiles")["query_target"]
    .agg(lambda values: ";".join(sorted(set(values))))
    .reset_index(name="active_targets")
)
active = active.merge(compound_targets, on="canonical_smiles", how="left")
active["has_v1"] = active["active_targets"].str.contains("TRPV1", regex=False)
active["has_non_v1"] = active["active_targets"].str.contains("TRPV2|TRPV3|TRPV4|TRPA1", regex=True)
selectivity = active[active["has_v1"] ^ active["has_non_v1"]].copy()
selectivity = selectivity.drop_duplicates("canonical_smiles")
selectivity["label"] = selectivity["has_v1"].astype(int)
selectivity.to_csv(V1_SELECTIVITY_TABLE, index=False)


########## 5. report ##########
print(f"Saved: {ACTIVITY_TABLE} ({len(activity)} rows)")
print(f"Saved: {TRPV1_BINDING_TABLE} ({len(trpv1_binding)} rows)")
print(f"Saved: {V1_SELECTIVITY_TABLE} ({len(selectivity)} rows)")

