from pathlib import Path
import time
import pandas as pd
import requests


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
TARGET_FILE = PROJECT_DIR / "configs/trpv_targets.tsv"
OUT_FILE = PROJECT_DIR / "outputs/ligands/chembl_trp_activity_raw.csv"
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
STANDARD_TYPES = {"IC50", "Ki", "Kd", "EC50"}
PAGE_LIMIT = 1000


########## 1. helpers ##########
def get_json(url, params=None):
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def find_target_ids(uniprot):
    data = get_json(
        f"{BASE_URL}/target_component.json",
        {"accession": uniprot, "limit": PAGE_LIMIT},
    )
    target_ids = []
    for component in data.get("target_components", []):
        for target in component.get("targets", []):
            target_id = target.get("target_chembl_id")
            if target_id:
                target_ids.append(target_id)
    return sorted(set(target_ids))


def fetch_activities(target_id):
    rows = []
    offset = 0
    while True:
        data = get_json(
            f"{BASE_URL}/activity.json",
            {
                "target_chembl_id": target_id,
                "limit": PAGE_LIMIT,
                "offset": offset,
            },
        )
        batch = data.get("activities", [])
        rows.extend(batch)
        if not data.get("page_meta", {}).get("next"):
            break
        offset += PAGE_LIMIT
        time.sleep(0.2)
    return rows


def classify_binder(row):
    pchembl = pd.to_numeric(row.get("pchembl_value"), errors="coerce")
    value = pd.to_numeric(row.get("standard_value"), errors="coerce")
    units = str(row.get("standard_units") or "").lower()
    if pd.notna(pchembl):
        if pchembl >= 6.0:
            return 1
        if pchembl <= 5.0:
            return 0
    if units == "nm" and pd.notna(value):
        if value <= 1000:
            return 1
        if value >= 10000:
            return 0
    return pd.NA


########## 2. fetch ##########
targets = pd.read_csv(TARGET_FILE, sep="\t")
all_rows = []
for _, target in targets.iterrows():
    target_ids = find_target_ids(target["uniprot"])
    for target_id in target_ids:
        for row in fetch_activities(target_id):
            if row.get("standard_type") not in STANDARD_TYPES:
                continue
            row["query_target"] = target["target"]
            row["query_uniprot"] = target["uniprot"]
            row["target_role"] = target["role"]
            all_rows.append(row)


########## 3. normalize ##########
df = pd.DataFrame(all_rows)
if df.empty:
    raise RuntimeError("No ChEMBL activity records were collected.")

keep_cols = [
    "query_target",
    "query_uniprot",
    "target_role",
    "target_chembl_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "assay_chembl_id",
    "assay_type",
    "document_chembl_id",
]
for col in keep_cols:
    if col not in df.columns:
        df[col] = pd.NA

df = df[keep_cols].drop_duplicates()
df["binder_label"] = df.apply(classify_binder, axis=1)
df["v1_label"] = (df["query_target"] == "TRPV1").astype(int)


########## 4. save ##########
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_FILE, index=False)
print(f"Saved: {OUT_FILE} ({len(df)} rows)")

