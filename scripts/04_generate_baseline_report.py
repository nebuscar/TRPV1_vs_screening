from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import load
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    PrecisionRecallDisplay,
    RocCurveDisplay,
    average_precision_score,
    confusion_matrix,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
MODEL_DIR = PROJECT_DIR / "outputs/models/ligand_baseline"
REPORT_DIR = PROJECT_DIR / "outputs/reports/ligand_baseline"
TASKS = {
    "trpv1_binding": PROJECT_DIR / "outputs/ligands/trpv1_binding_binary_table.csv",
    "v1_selectivity": PROJECT_DIR / "outputs/ligands/trpv1_vs_nonv1_selectivity_table.csv",
}
ACTIVITY_FILE = PROJECT_DIR / "outputs/ligands/trp_activity_model_table.csv"
RADIUS = 2
N_BITS = 2048
SEED = 42
TOP_N = 50


########## 1. helpers ##########
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
    return np.vstack(features), df.iloc[keep]["label"].astype(int).to_numpy(), df.iloc[keep].reset_index(drop=True)


def scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def score_smiles(smiles, model):
    rows = []
    valid = []
    for idx, item in enumerate(smiles):
        fp = morgan_fp(item)
        if fp is None:
            continue
        rows.append(fp)
        valid.append(idx)
    probs = model.predict_proba(np.vstack(rows))[:, 1]
    out = pd.Series(np.nan, index=range(len(smiles)), dtype=float)
    out.iloc[valid] = probs
    return out


########## 2. evaluation plots ##########
REPORT_DIR.mkdir(parents=True, exist_ok=True)
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
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for model_name, model in models.items():
        model.fit(x_train, y_train)
        prob = model.predict_proba(x_test)[:, 1]
        pred = (prob >= 0.5).astype(int)
        RocCurveDisplay.from_predictions(y_test, prob, name=model_name, ax=axes[0])
        PrecisionRecallDisplay.from_predictions(y_test, prob, name=model_name, ax=axes[1])
        metrics.append(
            {
                "task": task_name,
                "model": model_name,
                "roc_auc": float(roc_auc_score(y_test, prob)),
                "pr_auc": float(average_precision_score(y_test, prob)),
                "mcc": float(matthews_corrcoef(y_test, pred)),
            }
        )
    rf = models["random_forest"]
    rf_prob = rf.predict_proba(x_test)[:, 1]
    rf_pred = (rf_prob >= 0.5).astype(int)
    ConfusionMatrixDisplay(confusion_matrix(y_test, rf_pred)).plot(ax=axes[2], colorbar=False)
    axes[0].set_title(f"{task_name} ROC")
    axes[1].set_title(f"{task_name} PR")
    axes[2].set_title(f"{task_name} RF confusion")
    fig.tight_layout()
    fig.savefig(REPORT_DIR / f"{task_name}_evaluation.png", dpi=300)
    plt.close(fig)


########## 3. data summary plots ##########
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for ax, (task_name, input_file) in zip(axes, TASKS.items()):
    df = pd.read_csv(input_file)
    df["scaffold"] = df["canonical_smiles"].map(scaffold)
    df["label"].value_counts().sort_index().plot(kind="bar", ax=ax, color=["#4f6d7a", "#c1666b"])
    ax.set_title(f"{task_name} labels")
    ax.set_xlabel("label")
    ax.set_ylabel("count")
fig.tight_layout()
fig.savefig(REPORT_DIR / "label_distribution.png", dpi=300)
plt.close(fig)


########## 4. top candidates ##########
activity = pd.read_csv(ACTIVITY_FILE)
activity = activity.drop_duplicates("canonical_smiles").copy()
binding_model = load(MODEL_DIR / "trpv1_binding_random_forest.joblib")
selectivity_model = load(MODEL_DIR / "v1_selectivity_random_forest.joblib")
activity["trpv1_binding_score"] = score_smiles(activity["canonical_smiles"].tolist(), binding_model)
activity["v1_selectivity_score"] = score_smiles(activity["canonical_smiles"].tolist(), selectivity_model)
activity["combined_score"] = activity[["trpv1_binding_score", "v1_selectivity_score"]].mean(axis=1)
activity["scaffold"] = activity["canonical_smiles"].map(scaffold)
top = (
    activity.sort_values(["combined_score", "max_pchembl"], ascending=[False, False])
    .head(TOP_N)
    .reset_index(drop=True)
)
top.to_csv(REPORT_DIR / "top50_known_trpv1_candidate_compounds.csv", index=False)


########## 5. molecule grid ##########
mols = [Chem.MolFromSmiles(s) for s in top.head(20)["canonical_smiles"]]
legends = [
    f"{i+1}: B={b:.2f}, S={s:.2f}"
    for i, (b, s) in enumerate(zip(top.head(20)["trpv1_binding_score"], top.head(20)["v1_selectivity_score"]))
]
img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(260, 180), legends=legends)
img.save(str(REPORT_DIR / "top20_known_candidate_structures.png"))


########## 6. save ##########
(REPORT_DIR / "baseline_report_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
print(f"Saved: {REPORT_DIR}")

