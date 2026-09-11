from pathlib import Path
import os
import time
import numpy as np
import pandas as pd
import requests
from joblib import load
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Lipinski, QED, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
ACTIVITY_FILE = PROJECT_DIR / "outputs/ligands/trp_activity_model_table.csv"
MODEL_DIR = PROJECT_DIR / "outputs/models/ligand_baseline"
OUT_DIR = PROJECT_DIR / "outputs/reports/unknown_chembl_screen"
RAW_LIBRARY = PROJECT_DIR / "outputs/screening_library/chembl_druglike_unknown_raw.csv"
SCORED_LIBRARY = OUT_DIR / "chembl_druglike_unknown_scored.csv"
TOP_HITS = OUT_DIR / "top100_putative_new_hits.csv"
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PAGE_LIMIT = 1000
MAX_MOLECULES = int(os.environ.get("MAX_MOLECULES", "100000"))
RADIUS = 2
N_BITS = 2048


########## 1. helpers ##########
def get_json(url, params=None):
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()
    return response.json()


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def scaffold(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def morgan_fp(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def mol_props(smiles, catalog):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    qed = QED.qed(mol)
    pains = len(catalog.GetMatches(mol))
    lipinski_violations = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    return {
        "MW": mw,
        "LogP": logp,
        "HBD": hbd,
        "HBA": hba,
        "TPSA": tpsa,
        "RotBonds": rotb,
        "QED": qed,
        "PAINS_count": pains,
        "Lipinski_violations": lipinski_violations,
        "Veber_pass": rotb <= 10 and tpsa <= 140,
    }


def score_frame(df, binding_model, selectivity_model):
    rows = []
    valid = []
    for idx, smiles in enumerate(df["canonical_smiles"]):
        fp = morgan_fp(smiles)
        if fp is None:
            continue
        rows.append(fp)
        valid.append(idx)
    x = np.vstack(rows)
    scored = df.iloc[valid].copy()
    scored["trpv1_binding_score"] = binding_model.predict_proba(x)[:, 1]
    scored["v1_selectivity_score"] = selectivity_model.predict_proba(x)[:, 1]
    return scored


########## 2. fetch unknown library ##########
OUT_DIR.mkdir(parents=True, exist_ok=True)
RAW_LIBRARY.parent.mkdir(parents=True, exist_ok=True)
known = pd.read_csv(ACTIVITY_FILE)
known_smiles = set(known["canonical_smiles"].dropna().map(canonicalize).dropna())
known_scaffolds = set(known["canonical_smiles"].dropna().map(scaffold).dropna())

rows = []
offset = 0
while len(rows) < MAX_MOLECULES:
    data = get_json(
        BASE_URL,
        {
            "limit": PAGE_LIMIT,
            "offset": offset,
            "molecule_properties__full_mwt__lte": 500,
            "molecule_properties__alogp__lte": 5,
            "molecule_properties__hba__lte": 10,
            "molecule_properties__hbd__lte": 5,
        },
    )
    batch = data.get("molecules", [])
    if not batch:
        break
    for item in batch:
        structures = item.get("molecule_structures") or {}
        props = item.get("molecule_properties") or {}
        smiles = canonicalize(structures.get("canonical_smiles"))
        if not smiles or smiles in known_smiles:
            continue
        rows.append(
            {
                "molecule_chembl_id": item.get("molecule_chembl_id"),
                "canonical_smiles": smiles,
                "pref_name": item.get("pref_name"),
                "max_phase": item.get("max_phase"),
                "chembl_full_mwt": props.get("full_mwt"),
                "chembl_alogp": props.get("alogp"),
            }
        )
        if len(rows) >= MAX_MOLECULES:
            break
    if not data.get("page_meta", {}).get("next"):
        break
    offset += PAGE_LIMIT
    time.sleep(0.1)

library = pd.DataFrame(rows).drop_duplicates("canonical_smiles")
library.to_csv(RAW_LIBRARY, index=False)


########## 3. score ##########
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog(params)

binding_model = load(MODEL_DIR / "trpv1_binding_random_forest.joblib")
selectivity_model = load(MODEL_DIR / "v1_selectivity_random_forest.joblib")
scored = score_frame(library, binding_model, selectivity_model)
props = pd.DataFrame([mol_props(s, catalog) for s in scored["canonical_smiles"]])
scored = pd.concat([scored.reset_index(drop=True), props], axis=1)
scored["scaffold"] = scored["canonical_smiles"].map(scaffold)
scored["novel_scaffold_to_trp_train"] = ~scored["scaffold"].isin(known_scaffolds)
scored["combined_score"] = (
    0.50 * scored["trpv1_binding_score"]
    + 0.35 * scored["v1_selectivity_score"]
    + 0.15 * scored["QED"]
)
scored.to_csv(SCORED_LIBRARY, index=False)


########## 4. select top hits ##########
hits = scored[
    (scored["trpv1_binding_score"] >= 0.70)
    & (scored["v1_selectivity_score"] >= 0.70)
    & (scored["QED"] >= 0.35)
    & (scored["Lipinski_violations"] <= 1)
    & (scored["Veber_pass"])
    & (scored["PAINS_count"] == 0)
].copy()
hits = hits.sort_values(
    ["novel_scaffold_to_trp_train", "combined_score", "trpv1_binding_score", "v1_selectivity_score"],
    ascending=[False, False, False, False],
).head(100)
hits.to_csv(TOP_HITS, index=False)


########## 5. molecule grid ##########
top20 = hits.head(20)
mols = [Chem.MolFromSmiles(s) for s in top20["canonical_smiles"]]
legends = [
    f"{chembl_id}\\nB={b:.2f} S={s:.2f} QED={q:.2f}"
    for chembl_id, b, s, q in zip(
        top20["molecule_chembl_id"],
        top20["trpv1_binding_score"],
        top20["v1_selectivity_score"],
        top20["QED"],
    )
]
img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(320, 220), legends=legends)
img.save(str(OUT_DIR / "top20_putative_new_hit_structures.png"))


########## 6. report ##########
print(f"Saved: {RAW_LIBRARY} ({len(library)} rows)")
print(f"Saved: {SCORED_LIBRARY} ({len(scored)} rows)")
print(f"Saved: {TOP_HITS} ({len(hits)} rows)")

