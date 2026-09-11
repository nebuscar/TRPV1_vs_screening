from pathlib import Path
import subprocess
import numpy as np
import pandas as pd


########## 0. params ##########
PROJECT_DIR = Path("/home/nizhu/Projects/06.TRPV1_vs_screening")
RAW_DIR = PROJECT_DIR / "outputs/receptors"
OUT_DIR = PROJECT_DIR / "outputs/receptors"
DOCKING_PDB_DIR = OUT_DIR / "docking_pdb"
DOCKING_PDBQT_DIR = OUT_DIR / "docking_pdbqt"
MD_DIR = OUT_DIR / "md"
META_DIR = OUT_DIR / "metadata"
KNOWN_SITES = {
    "8X94": {"ligand": "EZI", "chain": "A", "site": "antagonist_vanilloid"},
    "7LPE": {"ligand": "4DY", "chain": "A", "site": "capsaicin_vanilloid"},
    "5IS0": {"ligand": "6ET", "chain": "B", "site": "capsazepine_vanilloid"},
}
PDB_IDS = ["11CJ", "5IS0", "6MHO", "6OT2", "7LPE", "8X94"]
WATER_NAMES = {"HOH", "WAT"}
ION_NAMES = {"NA", "CL", "K", "CA", "MG", "ZN"}
POCKET_RADIUS = 18.0
BOX_PADDING = 10.0


########## 1. helpers ##########
def parse_xyz(line):
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def residue_key(line):
    return line[21], line[22:26], line[26]


def fix_residue_name(line):
    if line.startswith("ATOM") and line[17:20] == "HIS":
        return line[:17] + "HIE" + line[20:]
    return line


def write_pdb(path, lines):
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def run_receptor_prep(input_pdb, output_base, center=None, size=None, allow_bad_res=False):
    command = [
        "mk_prepare_receptor.py",
        "--read_pdb",
        str(input_pdb),
        "--output_basename",
        str(output_base),
        "--write_pdbqt",
    ]
    if allow_bad_res:
        command += ["--delete_bad_res"]
    if center is not None and size is not None:
        command += [
            "--box_center",
            *(f"{v:.3f}" for v in center),
            "--box_size",
            *(f"{v:.3f}" for v in size),
        ]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, result.stderr[-1000:], "delete_bad_res" if allow_bad_res else "strict"


def run_receptor_prep_with_fallback(input_pdb, output_base, center=None, size=None):
    code, err, mode = run_receptor_prep(input_pdb, output_base, center, size)
    if code == 0:
        return code, err, mode
    return run_receptor_prep(input_pdb, output_base, center, size, allow_bad_res=True)


########## 2. prepare directories ##########
for folder in [DOCKING_PDB_DIR, DOCKING_PDBQT_DIR, MD_DIR, META_DIR]:
    folder.mkdir(parents=True, exist_ok=True)


########## 3. process receptors ##########
summary = []
for pdb_id in PDB_IDS:
    pdb_file = RAW_DIR / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        summary.append({"pdb_id": pdb_id, "raw_pdb": str(pdb_file), "status": "missing"})
        continue
    atom_lines = []
    md_lines = []
    ligand_lines = []
    n_water = 0
    n_ion = 0
    n_other_het = 0
    site = KNOWN_SITES.get(pdb_id)
    for line in pdb_file.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            fixed = fix_residue_name(line)
            atom_lines.append(fixed)
            md_lines.append(fixed)
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in WATER_NAMES:
                n_water += 1
                continue
            if resname in ION_NAMES:
                n_ion += 1
                md_lines.append(line)
                continue
            n_other_het += 1
            md_lines.append(line)
            if site and resname == site["ligand"] and line[21].strip() == site["chain"]:
                ligand_lines.append(line)
    protein_pdb = DOCKING_PDB_DIR / f"{pdb_id}_protein_only.pdb"
    md_pdb = MD_DIR / f"{pdb_id}_md_input_keep_nonwater_hetero.pdb"
    write_pdb(protein_pdb, atom_lines)
    write_pdb(md_pdb, md_lines)
    full_base = DOCKING_PDBQT_DIR / f"{pdb_id}_protein_only"
    full_code, full_err, full_mode = run_receptor_prep_with_fallback(protein_pdb, full_base)
    row = {
        "pdb_id": pdb_id,
        "raw_pdb": str(pdb_file),
        "protein_only_pdb": str(protein_pdb),
        "md_input_pdb": str(md_pdb),
        "n_protein_atoms": len(atom_lines),
        "n_water_removed": n_water,
        "n_ions_kept_for_md": n_ion,
        "n_nonwater_hetero_kept_for_md": n_other_het,
        "protein_pdbqt_status": full_code,
        "protein_pdbqt_mode": full_mode,
        "protein_pdbqt_error_tail": full_err,
    }
    if site and ligand_lines:
        ligand_xyz = np.array([parse_xyz(line) for line in ligand_lines])
        center = ligand_xyz.mean(axis=0)
        size = np.maximum(ligand_xyz.max(axis=0) - ligand_xyz.min(axis=0) + BOX_PADDING, 18.0)
        selected = {
            residue_key(line)
            for line in atom_lines
            if np.linalg.norm(parse_xyz(line) - center) <= POCKET_RADIUS
        }
        pocket_lines = [line for line in atom_lines if residue_key(line) in selected]
        pocket_pdb = DOCKING_PDB_DIR / f"{pdb_id}_{site['site']}_pocket.pdb"
        pocket_base = DOCKING_PDBQT_DIR / f"{pdb_id}_{site['site']}_pocket"
        write_pdb(pocket_pdb, pocket_lines)
        pocket_code, pocket_err, pocket_mode = run_receptor_prep_with_fallback(pocket_pdb, pocket_base, center, size)
        row.update(
            {
                "site_name": site["site"],
                "site_ligand": site["ligand"],
                "site_chain": site["chain"],
                "pocket_pdb": str(pocket_pdb),
                "pocket_pdbqt": str(pocket_base) + ".pdbqt",
                "pocket_pdbqt_status": pocket_code,
                "pocket_pdbqt_mode": pocket_mode,
                "pocket_pdbqt_error_tail": pocket_err,
                "box_center_x": center[0],
                "box_center_y": center[1],
                "box_center_z": center[2],
                "box_size_x": size[0],
                "box_size_y": size[1],
                "box_size_z": size[2],
                "n_pocket_atoms": len(pocket_lines),
                "n_site_ligand_atoms": len(ligand_lines),
            }
        )
    summary.append(row)


########## 4. save ##########
summary_df = pd.DataFrame(summary)
summary_file = META_DIR / "receptor_preprocessing_summary.csv"
summary_df.to_csv(summary_file, index=False)
print(f"Saved: {summary_file} ({len(summary_df)} rows)")

