# Deep Audit Report: Tubulin Virtual Screening Pipeline

> **Audit Date:** 2026-05-29
> **Auditor:** Antigravity Deep Analysis Engine
> **Scope:** Every source file, data file, result, published report, and documentation artifact in the project.
> **Verdict:** The core ML pipeline (`get_ic50.py`, `train_ml.py`, `resolve_and_predict.py`) is **scientifically sound**, but the project contains **several critical documentation mismatches**, an **ESOL formula bug**, a **legacy incompatible script**, and **synthetically fabricated MD simulation data**.

---

## Summary of Findings

| Severity | Count | Description |
|:---|:---:|:---|
| 🔴 CRITICAL | 3 | MD data fabrication, ESOL formula errors, README/report showing OLD pipeline data |
| 🟠 MAJOR | 4 | Legacy script incompatibility, missing docking scores, Table 2 Pa/Pi values hardcoded, ChEMBL target ambiguity |
| 🟡 MINOR | 3 | Feature selection leakage (minor), control SMILES not PubChem-canonical, solubility unit conversion |
| 🟢 VERIFIED CORRECT | 12 | PubChem CIDs, Morgan FP generation, MinMaxScaler usage, Tanimoto penalty, Lipinski calculations, etc. |

---

## 🔴 CRITICAL ISSUES

### CRITICAL-1: MD Trajectory Data is FABRICATED (Not Real Simulation Output)

**File:** [analyze_md_trajectory.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/analyze_md_trajectory.py)

The function `simulate_md_data()` (line 11) generates **entirely synthetic/fake** trajectory data using NumPy random number generators. It does NOT read any actual GROMACS `.xvg`, `.trr`, or `.xtc` trajectory files.

```python
# Line 20: backbone RMSD is fabricated with exponential decay + noise
backbone_rmsd = 0.5 + 1.7 * (1.0 - np.exp(-time_ns / 8.0)) + np.random.normal(0, 0.05, len(time_ns))

# Line 24: ligand RMSD is also fabricated
withaferin_rmsd = 0.3 + 1.2 * (1.0 - np.exp(-time_ns / 5.0)) + np.random.normal(0, 0.06, len(time_ns))
```

The GROMACS shell script [run_md_gromacs.sh](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/run_md_gromacs.sh) is a **template** that was never executed. No actual MD simulation was performed.

> [!CAUTION]
> **Scientific Integrity Risk:** The MD trajectory plots (`md_rmsd_rmsf.png`) and CSV data (`md_simulation_results.csv`, `md_rmsf_results.csv`) are presented as simulation results but are **mathematically generated** data. If used in a publication or thesis, this would constitute data fabrication. The plots should be explicitly labeled as "Illustrative/Theoretical Framework" or the simulation should be actually run.

**Impact:** The dashboard's MD Simulation tab and the published web page both present this fabricated data as real results.

---

### CRITICAL-2: ESOL (Delaney Solubility) Formula Has TWO Coefficient Errors

**File:** [admet_prediction.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/admet_prediction.py) — Line 34

**The correct Delaney ESOL equation (from the 2004 paper):**
```
LogS = 0.16 - 0.63*cLogP - 0.0062*MW + 0.066*RB - 0.74*AP
```

**What the code actually implements:**
```python
logs = 0.16 - (0.63 * logp) - (0.0062 * mw) - (0.0062 * rtb) + (0.05 * aromatic_prop)
```

| Parameter | Correct Coefficient | Coded Coefficient | Error |
|:---|:---:|:---:|:---|
| Rotatable Bonds (RB) | **+0.066** | **-0.0062** | Wrong sign AND wrong value (should be +0.066, not -0.0062) |
| Aromatic Proportion (AP) | **-0.74** | **+0.05** | Wrong sign AND wrong magnitude (off by 14.8×) |

> [!WARNING]
> This means ALL estimated solubility values in `admet_predictions.csv` and the ADMET table in the dashboard/README are **incorrect**. The error particularly affects compounds with many rotatable bonds or high aromatic content.

---

### CRITICAL-3: README.md and drug_discovery_report.md Show OLD/Stale Pipeline Data

**Files:**
- [README.md](file:///C:/Users/k4ran/OneDrive/Desktop/ML/README.md)
- [drug_discovery_report.md](file:///C:/Users/k4ran/OneDrive/Desktop/ML/results/drug_discovery_report.md)

Both documents describe the **original, pre-overhaul pipeline** and were never updated to reflect the current code:

| Claim in README/Report | Actual Current Code | Status |
|:---|:---|:---:|
| "60 compounds (30 active, 30 inactive)" | **300 compounds (150 active, 150 inactive)** via `get_ic50.py` | ❌ WRONG |
| "57 unique compounds" | **300 unique compounds** | ❌ WRONG |
| "MACCS Keys (166 bits)" | **Morgan ECFP4 (2048 bits)** via `AllChem.GetMorganFingerprintAsBitVect` | ❌ WRONG |
| "383 dimensions" feature matrix | **2048 + ~200 RDKit descriptors = ~2248**, then SelectKBest(k=100) → **100 features** | ❌ WRONG |
| "StandardScaler" | **MinMaxScaler** (line 167 of `train_ml.py`) | ❌ WRONG |
| Performance Table (RF AUC=0.8456, LR AUC=0.8444) | **RF AUC=0.8144, LR AUC=0.7880, SVM AUC=0.8204, XGBoost AUC=0.8269** (from `exact_tables_report.md`) | ❌ WRONG |
| Hit Predictions: "All 7 predicted Inactive" (consensus prob 0.26-0.42) | **5 predicted Active** (consensus prob 0.50-0.57) in current `hit_predictions_results.csv` | ❌ WRONG |
| Feature importances: "MACCS_66, MACCS_158, MACCS_112" | Current model uses **ECFP4 bits**, not MACCS keys | ❌ WRONG |

> [!CAUTION]
> The README is the public face of the GitHub repository. Anyone reading it will get a completely incorrect understanding of the pipeline's methodology, training data size, feature engineering approach, and results.

---

## 🟠 MAJOR ISSUES

### MAJOR-1: Legacy `predict_hits.py` is Incompatible with Current Trained Model

**File:** [predict_hits.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/predict_hits.py)

This script uses **MACCS keys** (line 28: `MACCSkeys.GenMACCSKeys`) and does NOT apply:
- Morgan fingerprints
- SelectKBest feature selection
- The 0.05 physical descriptor downweighting
- Tanimoto Vinca scaffold penalty

If someone runs this script with the current `trained_models.pkl` (which was trained on Morgan FP + SelectKBest=100 features), it would crash or produce completely garbage predictions because the feature dimensions don't match.

**The correct prediction script is** [resolve_and_predict.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/resolve_and_predict.py), which properly uses Morgan FPs, MinMaxScaler, selector, and Vinca penalty.

**Recommendation:** Delete or archive `predict_hits.py` to avoid confusion.

---

### MAJOR-2: Docking Scores Are Not Used in the ML Prediction Pipeline

The user's AutoDock Vina docking scores (Withanolide A = -10.6, Withaferin A = -10.2, etc.) are **not used as features** in the ML classification. The pipeline only uses:
1. Morgan Fingerprints (ECFP4)
2. RDKit 2D physical descriptors (scaled and downweighted)

The docking scores are only mentioned in the published documentation and dashboard as context, but they do NOT influence the ML model's Active/Inactive predictions. This is **not a bug** — it's actually correct methodology (docking scores are pre-screening, ML is post-screening) — but it should be clearly documented.

---

### MAJOR-3: Table 2 Pa/Pi Values Are HARDCODED, Not Computed

**File:** [generate_tables.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/generate_tables.py) — Lines 140-181

The biological activity predictions (Pa for "Tubulin antagonist", "Apoptosis agonist", "Cytostatic", etc.) shown in Table 2 of `exact_tables_report.md` are **hardcoded dictionaries**, not computed from any model or PASS server:

```python
class_activities = {
    "Withaferin A": [
        ("Tubulin antagonist", 0.425, 0.025),  # <-- hardcoded
        ("Apoptosis agonist", 0.650, 0.012),    # <-- hardcoded
        ...
    ],
```

These Pa/Pi values appear to be **manually curated from literature or the PASS online server**, but there is no attribution or documentation of where exactly they came from. They are presented in the README as "PASS Biological Activity Predictions" but no PASS computation was performed in the pipeline.

> [!IMPORTANT]
> If these values are from an external PASS server run, that should be cited. If they are approximate literature estimates, they should be labeled as such.

---

### MAJOR-4: ChEMBL Target ID Ambiguity

**File:** [get_ic50.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/get_ic50.py) — Lines 54-58

The code queries **16 different ChEMBL targets**, but:
- **CHEMBL3394** = Tubulin beta chain (*Bos taurus*, bovine, SINGLE PROTEIN) — NOT human
- **CHEMBL3885647** = Tubulin alpha-1A/beta chains (PROTEIN COMPLEX GROUP) — general tubulin, not specifically beta-III
- **CHEMBL2597** = Tubulin beta-3 chain (TUBB3) — THIS is the actual human beta-III tubulin target

The query mixes bovine tubulin, general tubulin complexes, and multiple different tubulin subtypes. While this is acceptable for building a diverse training set (tubulin is highly conserved across species), the README/report should clarify that the training data comes from **multiple tubulin subtypes across species**, not exclusively "human β-III tubulin".

---

## 🟡 MINOR ISSUES

### MINOR-1: Feature Selection Fitted on Full Dataset Before CV (Minor Leakage)

**File:** [train_ml.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/train_ml.py) — Lines 178-179

```python
global_selector = SelectKBest(score_func=f_classif, k=100)
X_selected = global_selector.fit_transform(X, y)
```

The `global_selector` is fit on the **entire dataset** (all 300 compounds). This is then used for control validation (line 258). However, the cross-validation loop correctly uses **fold-level** feature selection (lines 211-213):

```python
fold_selector = SelectKBest(score_func=f_classif, k=100)
X_train_sel = fold_selector.fit_transform(X_train, y_train)
X_val_sel = fold_selector.transform(X_val)
```

So the **CV metrics are correct** (no leakage there). The global selector is only used for the final model and control validation, which is acceptable practice for production deployment.

**Verdict:** Minor issue, not a real problem. CV is done correctly.

---

### MINOR-2: Solubility mg/L Conversion Formula

**File:** [admet_prediction.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/admet_prediction.py) — Line 111

```python
solubility_mg_l = (10 ** logs) * mw * 1000.0  # convert mol/L to mg/L
```

ESOL LogS is in log(mol/L). So `10**logs` gives mol/L. To convert to mg/L: `mol/L × MW(g/mol) × 1000(mg/g)`. This gives mg/L, which is correct. However, since the LogS values themselves are wrong (CRITICAL-2), these solubility values are also wrong.

---

### MINOR-3: MinMaxScaler Fitted Before Train/Test Split

**File:** [train_ml.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/train_ml.py) — Lines 167-168

```python
scaler = MinMaxScaler()
rdkit_feats_scaled = scaler.fit_transform(rdkit_feats)
```

The scaler is fitted on the **entire dataset** before the CV split. Technically, for strict ML best practice, the scaler should be fitted only on the training fold. However:
1. The 0.05 downweighting makes physical features near-zero anyway
2. The feature selection further filters features
3. This is standard practice for small molecular datasets in cheminformatics

**Verdict:** Minor, acceptable.

---

## 🟢 VERIFIED CORRECT

### ✅ VERIFIED-1: PubChem Compound IDs
All 7 CIDs verified against PubChem:

| Compound | Claimed CID | Verified | Status |
|:---|:---:|:---:|:---:|
| Withaferin A | 265237 | 265237 | ✅ CORRECT |
| 27-O-acetyl-Withaferin A | 57328756 | 57328756 | ✅ CORRECT |
| Withanolide A | 11294368 | 11294368 | ✅ CORRECT |
| Celastrol | 122724 | 122724 | ✅ CORRECT |
| Pristimerin | 159516 | 159516 | ✅ CORRECT |
| Tingenone | 101520 | 101520 | ✅ CORRECT |
| Alpha-Glycyrrhizin | 158471 | 158471 | ✅ CORRECT |

### ✅ VERIFIED-2: Morgan Fingerprint Generation
`train_ml.py` line 39: `AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)` — This correctly generates ECFP4 fingerprints (Morgan radius=2 is equivalent to ECFP diameter=4).

### ✅ VERIFIED-3: IC50 Threshold Classification
- Active: IC50 < 1 µM (1000 nM) — standard for potent inhibitors
- Inactive: IC50 > 10 µM (10000 nM) — standard for weak/inactive
- Gap between 1-10 µM is excluded (ambiguous zone) — correct practice

### ✅ VERIFIED-4: Lipinski Rule of Five Implementation
`admet_prediction.py` lines 103-107:
- MW > 500 ✅ (standard)
- LogP > 5 ✅ (standard)
- HBD > 5 ✅ (standard)
- HBA > 10 ✅ (standard)

### ✅ VERIFIED-5: Tanimoto Vinca Scaffold Penalty
The penalty `prob * (1 - max_vinca_sim)` applied when similarity > 0.70 is scientifically justified. Vincristine has high Tanimoto similarity to Vinblastine/Vinorelbine/Vindesine (all Vinca alkaloids), and penalizing their active probability at the taxane site is correct because they are mechanistically inactive at the taxane pocket.

### ✅ VERIFIED-6: Cross-Validation Is Done Correctly
`train_ml.py` lines 206-226: Fold-level SelectKBest is applied independently to each fold, preventing feature selection leakage in the CV performance estimates.

### ✅ VERIFIED-7: Control Drug Validation Logic
Paclitaxel (taxane stabilizer → Active) and Vincristine (Vinca destabilizer → Inactive at taxane site) are correctly excluded from training (line 12: `EXCLUDE_CHEMBL_IDS`) and used as external validation controls.

### ✅ VERIFIED-8: Model Hyperparameters Are Reasonable
- SVM C=0.3 (moderate regularization) ✅
- RF max_depth=3 (shallow trees prevent overfitting on small datasets) ✅
- LR C=0.005 (strong L2 regularization) ✅
- XGBoost max_depth=2, reg_alpha=1.0, reg_lambda=3.0 (heavily regularized) ✅

### ✅ VERIFIED-9: Paclitaxel/Vincristine SMILES Are Valid
Both control drug SMILES parse correctly in RDKit and match the known structures.

### ✅ VERIFIED-10: Pharmacological Rationale for Inactive Predictions
The explanation that Withaferin A, Celastrol, and Pristimerin act via covalent cysteine modification (Cys239, Cys303) rather than taxane pocket binding is **confirmed by literature** (NIH/PubMed sources verified).

### ✅ VERIFIED-11: Feature Importance Extraction
`feature_importance.py` correctly extracts RF `feature_importances_` and LR `coef_[0]` from the trained models. Feature name mapping using the SelectKBest support mask is correct.

### ✅ VERIFIED-12: PAINS Filter
`admet_prediction.py` lines 37-47: Correctly uses RDKit's `FilterCatalog` with `PAINS` filter catalog. Implementation is standard.

---

## 📋 RECOMMENDED FIXES (Priority Order)

### Priority 1: Fix ESOL Formula (CRITICAL-2)

In [admet_prediction.py](file:///C:/Users/k4ran/OneDrive/Desktop/ML/src/admet_prediction.py), line 34, change:
```diff
-    logs = 0.16 - (0.63 * logp) - (0.0062 * mw) - (0.0062 * rtb) + (0.05 * aromatic_prop)
+    logs = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rtb) - (0.74 * aromatic_prop)
```

### Priority 2: Update README.md (CRITICAL-3)

Rewrite the README to accurately reflect:
- 300-compound balanced dataset (150 active/150 inactive)
- Morgan ECFP4 fingerprints (2048 bits) + RDKit 2D descriptors
- MinMaxScaler with 0.05 downweighting
- SelectKBest(k=100) feature selection
- Current CV performance numbers (SVM AUC=0.82, RF AUC=0.81, LR AUC=0.79, XGBoost AUC=0.83)
- Current hit predictions (5 Active, 2 Inactive)

### Priority 3: Label MD Simulation as Illustrative (CRITICAL-1)

Either:
- **(a)** Clearly label all MD data/plots as "Illustrative/Theoretical Framework — No actual GROMACS simulation was performed"
- **(b)** Actually run the GROMACS simulation using `run_md_gromacs.sh`

### Priority 4: Remove/Archive Legacy predict_hits.py (MAJOR-1)

Delete or rename to `predict_hits_LEGACY.py` to prevent confusion.

### Priority 5: Document Pa/Pi Source (MAJOR-3)

Add a comment in `generate_tables.py` citing the PASS Online server or literature source for the biological activity predictions.

### Priority 6: Update drug_discovery_report.md (CRITICAL-3)

Same updates as README — align with current pipeline methodology and results.

---

## 🔍 External Database Verification Sources

| Database | Verification Target | Status |
|:---|:---|:---:|
| PubChem (NIH) | All 7 compound CIDs and structures | ✅ Verified |
| ChEMBL (EBI) | CHEMBL3394 = Bovine tubulin beta | ✅ Verified (but note: bovine, not human) |
| ChEMBL (EBI) | CHEMBL3885647 = Tubulin alpha/beta complex group | ✅ Verified (but note: not specifically beta-III) |
| ChEMBL (EBI) | CHEMBL2597 = Human TUBB3 (beta-III tubulin) | ✅ Verified |
| PubMed/NIH | Withaferin A covalent mechanism (Cys303/Cys239) | ✅ Verified |
| PubMed/NIH | Celastrol cysteine-reactive tubulin destabilizer | ✅ Verified |
| Delaney 2004 (ACS) | ESOL LogS formula coefficients | ✅ Verified (code has errors) |
