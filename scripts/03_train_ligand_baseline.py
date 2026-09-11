from pathlib import Path
import json
import numpy as np
import pandas as pd
from joblib import dump
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import train_test_split


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
TASKS = {
    "trpv1_binding": PROJECT_DIR / "outputs/ligands/trpv1_binding_binary_table.csv",
    "v1_selectivity": PROJECT_DIR / "outputs/ligands/trpv1_vs_nonv1_selectivity_table.csv",
}
OUT_DIR = PROJECT_DIR / "outputs/models/ligand_baseline"
RADIUS = 2
N_BITS = 2048
SEED = 42


########## 1. featurize ##########
def morgan_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=N_BITS)
    arr = np.zeros((N_BITS,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def build_matrix(df):
    features = []
    keep = []
    for idx, smiles in enumerate(df["canonical_smiles"]):
        fp = morgan_fp(smiles)
        if fp is None:
            continue
        features.append(fp)
        keep.append(idx)
    x = np.vstack(features)
    y = df.iloc[keep]["label"].astype(int).to_numpy()
    return x, y, df.iloc[keep].reset_index(drop=True)


########## 2. train ##########
def evaluate(name, model, x_test, y_test):
    prob = model.predict_proba(x_test)[:, 1]
    pred = (prob >= 0.5).astype(int)
    return {
        "task": name,
        "model": model.__class__.__name__,
        "roc_auc": float(roc_auc_score(y_test, prob)),
        "pr_auc": float(average_precision_score(y_test, prob)),
        "mcc": float(matthews_corrcoef(y_test, pred)),
    }


OUT_DIR.mkdir(parents=True, exist_ok=True)
metrics = []
for task_name, input_file in TASKS.items():
    df = pd.read_csv(input_file)
    x, y, used = build_matrix(df)
    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=0.2,
        random_state=SEED,
        stratify=y,
    )
    models = {
        "logistic_regression": LogisticRegression(max_iter=2000, class_weight="balanced"),
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            random_state=SEED,
            class_weight="balanced_subsample",
            n_jobs=-1,
        ),
    }
    for model_name, model in models.items():
        model.fit(x_train, y_train)
        metrics.append(evaluate(task_name, model, x_test, y_test))
        dump(model, OUT_DIR / f"{task_name}_{model_name}.joblib")
    used.to_csv(OUT_DIR / f"{task_name}_used_compounds.csv", index=False)


########## 3. save ##########
metrics_file = OUT_DIR / "metrics.json"
metrics_file.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
print(f"Saved: {metrics_file}")

