from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
PDB_ID = "8X94"
LIGAND_CODE = "EZI"
LIGAND_CHAIN = "A"
RAW_PDB = PROJECT_DIR / f"outputs/receptors/{PDB_ID}.pdb"
PREP_DIR = PROJECT_DIR / "outputs/receptors"
DOCK_DIR = PROJECT_DIR / f"outputs/docking/{PDB_ID.lower()}_top_hits"
TOP_HITS = PROJECT_DIR / "outputs/reports/unknown_chembl_screen/top100_putative_new_hits.csv"
PROTEIN_PDB = PREP_DIR / "docking_pdb/8X94_antagonist_vanilloid_pocket.pdb"
RECEPTOR_PREFIX = PREP_DIR / "docking_pdbqt/8X94_antagonist_vanilloid_pocket"
RECEPTOR_PDBQT = PREP_DIR / "docking_pdbqt/8X94_antagonist_vanilloid_pocket.pdbqt"
OUT_FILE = DOCK_DIR / "top20_vina_scores.csv"
POSE_DIR = DOCK_DIR / "poses"
N_DOCK = 20
BOX_PADDING = 10.0
POCKET_RADIUS = 18.0
EXHAUSTIVENESS = 8
N_POSES = 3
SEED = 42


########## 1. receptor and box ##########
PREP_DIR.mkdir(parents=True, exist_ok=True)
(PREP_DIR / "docking_pdb").mkdir(parents=True, exist_ok=True)
(PREP_DIR / "docking_pdbqt").mkdir(parents=True, exist_ok=True)
DOCK_DIR.mkdir(parents=True, exist_ok=True)
POSE_DIR.mkdir(parents=True, exist_ok=True)

ligand_xyz = []
atom_records = []
for line in RAW_PDB.read_text(errors="ignore").splitlines():
    if line.startswith("ATOM"):
        atom_records.append(line)
    if (
        line.startswith("HETATM")
        and line[17:20].strip() == LIGAND_CODE
        and line[21].strip() == LIGAND_CHAIN
    ):
        ligand_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

if not ligand_xyz:
    raise RuntimeError(f"No ligand {LIGAND_CODE} found in {RAW_PDB}")

xyz = np.array(ligand_xyz)
center = xyz.mean(axis=0)
size = np.maximum(xyz.max(axis=0) - xyz.min(axis=0) + BOX_PADDING, 18.0)

selected_residues = set()
for line in atom_records:
    coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if np.linalg.norm(coord - center) <= POCKET_RADIUS:
        selected_residues.add((line[21], line[22:26], line[26]))

protein_lines = []
for line in atom_records:
    residue_key = (line[21], line[22:26], line[26])
    if residue_key not in selected_residues:
        continue
    if line[17:20] == "HIS":
        line = line[:17] + "HIE" + line[20:]
    protein_lines.append(line)

PROTEIN_PDB.write_text("\n".join(protein_lines) + "\nEND\n")

if not RECEPTOR_PDBQT.exists():
    subprocess.run(
        [
            "mk_prepare_receptor.py",
            "--read_pdb",
            str(PROTEIN_PDB),
            "--output_basename",
            str(RECEPTOR_PREFIX),
            "--write_pdbqt",
            "--box_center",
            *(f"{v:.3f}" for v in center),
            "--box_size",
            *(f"{v:.3f}" for v in size),
        ],
        check=True,
    )


########## 2. ligand preparation ##########
def ligand_to_pdbqt(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, randomSeed=SEED, maxAttempts=1000)
    if status != 0:
        status = AllChem.EmbedMolecule(mol, randomSeed=SEED, maxAttempts=1000, useRandomCoords=True)
    if status != 0:
        raise ValueError("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    return PDBQTWriterLegacy.write_string(setups[0])[0]


########## 3. dock ##########
hits = pd.read_csv(TOP_HITS).head(N_DOCK).copy()
v = Vina(sf_name="vina", cpu=4, seed=SEED)
v.set_receptor(str(RECEPTOR_PDBQT))
v.compute_vina_maps(center=center.tolist(), box_size=size.tolist())

rows = []
for idx, row in hits.iterrows():
    chembl_id = row["molecule_chembl_id"]
    try:
        pdbqt = ligand_to_pdbqt(row["canonical_smiles"])
        v.set_ligand_from_string(pdbqt)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        pose_file = POSE_DIR / f"{idx + 1:03d}_{chembl_id}_vina.pdbqt"
        v.write_poses(str(pose_file), n_poses=N_POSES, overwrite=True)
        rows.append(
            {
                "rank_input": idx + 1,
                "molecule_chembl_id": chembl_id,
                "vina_best_kcal_mol": float(energies[0][0]),
                "vina_inter_kcal_mol": float(energies[0][1]),
                "trpv1_binding_score": row["trpv1_binding_score"],
                "v1_selectivity_score": row["v1_selectivity_score"],
                "QED": row["QED"],
                "MW": row["MW"],
                "LogP": row["LogP"],
                "pose_file": str(pose_file),
                "canonical_smiles": row["canonical_smiles"],
            }
        )
    except Exception as exc:
        rows.append(
            {
                "rank_input": idx + 1,
                "molecule_chembl_id": chembl_id,
                "vina_best_kcal_mol": np.nan,
                "vina_inter_kcal_mol": np.nan,
                "error": str(exc),
                "canonical_smiles": row["canonical_smiles"],
            }
        )


########## 4. save ##########
out = pd.DataFrame(rows).sort_values("vina_best_kcal_mol", na_position="last")
out.to_csv(OUT_FILE, index=False)
print(f"Saved: {OUT_FILE} ({len(out)} rows)")
print(f"Box center: {center.round(3).tolist()}")
print(f"Box size: {size.round(3).tolist()}")

