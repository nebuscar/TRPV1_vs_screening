# TRPV1 小分子虚拟筛选模型

本项目目标是围绕 TRPV1 建立可复现的小分子虚拟筛选流程，同时显式纳入 TRPV2、TRPV3、TRPV4、TRPA1 作为非 V1 同源/相关靶点对照。

## 建模目标

########## 0. target definition ##########

核心靶点为 TRPV1。非 V1 对照靶点包括 TRPV2、TRPV3、TRPV4、TRPA1。

模型需要解决两个问题：

1. 区分 V1 与非 V1 靶点相关配体/结合模式；
2. 对十万级化合物库进行批量筛选，综合蛋白动态构象特征、化合物结构特征和蛋白-配体结合模式，输出结合活性打分，并按 cutoff 判定结合或不结合。

## 推荐技术路线

########## 1. receptor data ##########

收集 TRPV1、TRPV2、TRPV3、TRPV4、TRPA1 的 PDB、AlphaFold、OPM/PPM 结构信息。TRPV 家族是膜蛋白，后续涉及动态构象时应优先使用膜环境 MD，而不是普通水盒体系。

########## 2. ligand activity data ##########

从 ChEMBL、BindingDB、PubChem BioAssay 获取各靶点小分子活性数据。优先保留标准化的 IC50、Ki、Kd、EC50，并统一到 nM 与 pActivity。

推荐初始标签：

- `binder = 1`：`pchembl_value >= 6.0` 或活性值 `<= 1000 nM`
- `binder = 0`：有明确 inactive/弱活性证据，或活性值 `>= 10000 nM`
- 中间区间作为 gray zone，训练二分类时先排除
- `target_class = V1` 或 `non_V1`

########## 3. baseline model ##########

先建立 ligand-only baseline：ECFP/Morgan fingerprint + RDKit descriptors + scaffold split。该模型回答“分子结构是否像 TRPV1 活性配体”，不声称识别蛋白动态构象。

推荐基线模型：

- Logistic Regression / Random Forest / XGBoost / LightGBM
- Chemprop D-MPNN
- scaffold split 外部验证
- 指标：ROC-AUC、PR-AUC、MCC、EF1%、EF5%、BEDROC

########## 4. structure-aware model ##########

对 TRPV1 与非 V1 结构做受体准备，围绕 vanilloid pocket 等口袋进行 docking 或 AI docking，提取蛋白-配体相互作用特征：

- docking score
- pose cluster
- hydrogen bond / hydrophobic / pi interaction
- pocket residue contact fingerprint
- ligand strain / torsion penalty
- TRPV1 vs 非 V1 的 score gap

推荐顺序：

1. AutoDock Vina/smina 建立可解释基线；
2. GNINA 或 DiffDock 增加深度学习打分；
3. 对 top hits 做 MD 稳定性验证。

########## 5. dynamic conformation features ##########

动态构象特征不建议一开始端到端硬学。更稳妥做法是先从 MD 轨迹中提取固定维度特征，再进入二分类/排序模型：

- pocket volume
- key residue distance
- RMSD/RMSF
- ligand contact occupancy
- water bridge occupancy
- binding pose stability
- MM/PBSA 或 MM/GBSA 近似能量

########## 6. final decision ##########

最终筛选不要只看一个分数。建议使用分层 gating：

1. 分子标准化与 PAINS/BRENK/无效结构过滤；
2. ligand-only TRPV1-like probability；
3. TRPV1 docking/AI docking score；
4. TRPV1 对非 V1 选择性 score gap；
5. ADMET 和 drug-likeness；
6. top hits 做 MD 稳定性复核。

## 目录约定

```text
data/
  raw/
    receptors/
    ligands/
    screening_library/
  processed/
    receptors/
    ligands/
    features/
  results/
    docking/
    models/
    reports/
configs/
scripts/
docs/
```

## 当前入口

- `scripts/00_initialize_project.py`：初始化目录；
- `scripts/01_fetch_chembl_trp_activity.py`：按 UniProt ID 从 ChEMBL 抓取 TRP 家族活性记录；
- `configs/trpv_targets.tsv`：TRP 靶点定义；
- `configs/screening_config.yaml`：初始 cutoff 与筛选配置。

## RCSB PDB 受体预处理

########## 7. receptor preprocessing ##########

RCSB PDB 原始结构不直接作为 docking 输入。当前项目已增加正式批量预处理脚本：

- `scripts/08_prepare_receptors.py`

该脚本从 `data/raw/receptors/*.pdb` 读取原始 PDB，输出两类受体文件：

- docking 输入：去除水分子和非目标 hetero，仅保留 protein，生成 `data/processed/receptors/docking_pdb/*_protein_only.pdb` 与 `data/processed/receptors/docking_pdbqt/*_protein_only.pdbqt`；
- known-ligand pocket 输入：对 `8X94/EZI`、`7LPE/4DY`、`5IS0/6ET` 按共晶配体坐标截取 vanilloid pocket，并生成 pocket PDBQT；
- MD 初始输入：删除水，但保留 protein、非水 hetero 和离子，写入 `data/processed/receptors/md/*_md_input_keep_nonwater_hetero.pdb`，后续应再接 OPM/CHARMM-GUI 膜体系构建。

当前预处理 summary：

- summary 文件：`data/processed/receptors/metadata/receptor_preprocessing_summary.csv`
- 全蛋白 PDBQT：`11CJ`、`5IS0`、`6MHO`、`6OT2`、`7LPE`、`8X94` 均已生成；
- vanilloid pocket PDBQT：`8X94_antagonist_vanilloid_pocket`、`7LPE_capsaicin_vanilloid_pocket`、`5IS0_capsazepine_vanilloid_pocket` 均已生成；
- `8X94` 通过严格模式生成；`5IS0/7LPE/6MHO/6OT2` 含部分模板不匹配或不完整残基，脚本会先严格转换，失败后用 `--delete_bad_res` 容错重试，并在 summary 中记录实际模式。

因此，当前项目已经有 PDB 结构清洗/受体制备步骤；后续 docking 不应直接使用 `data/raw/receptors/*.pdb`，而应使用 `data/processed/receptors/docking_pdbqt/*.pdbqt`。
