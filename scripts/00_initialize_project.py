from pathlib import Path


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
DIRS = [
    "outputs/receptors",
    "outputs/receptors/docking_pdb",
    "outputs/receptors/docking_pdbqt",
    "outputs/receptors/md",
    "outputs/receptors/metadata",
    "outputs/ligands",
    "outputs/features",
    "outputs/screening_library",
    "outputs/docking",
    "outputs/models",
    "outputs/reports",
    "logs",
]


########## 1. initialize ##########
for rel_path in DIRS:
    (PROJECT_DIR / rel_path).mkdir(parents=True, exist_ok=True)


########## 2. report ##########
print(f"Initialized: {PROJECT_DIR}")

