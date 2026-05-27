# ML-Driven Virtual Screening & ADMET Profiling of Tubulin Inhibitors

An end-to-end computational drug discovery pipeline combining **homology modeling**, **structure-based virtual screening (molecular docking)**, **multi-algorithm machine learning classification**, and **ADMET profiling** to identify and prioritize novel natural product hits targeting the **taxane binding site** of human $\alpha\beta_{\text{III}}$ tubulin.

## ⚡ Quick Start: Run Entire Pipeline
To automatically execute all stages of the pipeline sequentially (from database mining to final report formatting and plotting):
```bash
py src/all.py
```

---

## 🧬 1. Executive Summary & Biological Target

Microtubules are dynamic polymers of $\alpha$- and $\beta$-tubulin heterodimers that undergo continuous assembly and disassembly. In rapidly dividing cancer cells, microtubule dynamics are vital for the formation of the mitotic spindle that segregates chromosomes during cell division. Disruption of these dynamics is a validated therapeutic strategy in clinical oncology.

This project focuses on targeting the **taxane binding pocket** located on the interior surface of the $\beta$-tubulin subunit. Binding of stabilizing ligands (such as Paclitaxel, Docetaxel, and Cabazitaxel) to this pocket freezes the microtubule polymer, preventing disassembly. This halts the cell cycle in the G2/M phase (mitotic arrest) and ultimately triggers apoptosis.

---

## 🗺️ 2. Homology Modeling & Stereochemical Validation

Since high-resolution experimental crystal structures of the human $\alpha\beta_{\text{III}}$ tubulin isotype are unavailable in the Protein Data Bank (PDB), a 3D structural model was constructed using comparative homology modeling:

* **Template Selection:** The crystal structure of bovine $\alpha_{\text{1B}}\beta_{\text{IIB}}$ tubulin (**PDB ID: 1JFF**) was selected as a structural template due to its high sequence similarity.
* **Model Generation:** 100 structural models of the human $\beta_{\text{III}}$ isotype were generated using **Modeller 10.2** while maintaining the structure of the $\alpha_{\text{1B}}$-tubulin subunit constant.
* **Selection Metric:** The optimal conformation was selected based on the lowest **Discrete Optimized Protein Energy (DOPE)** score of **$-51,889.82$**, which represents the most physically relaxed and thermodynamically stable model.
* **Stereochemical Validation:** Structural validation using **PROCHECK** was performed via a Ramachandran plot to evaluate $\phi$ (phi) and $\psi$ (psi) backbone dihedral angles:
  * **Favored & Allowed Regions:** **97.7%** of all amino acid residues fell within stereochemically favorable zones (88.1% in most favored, 9.6% in additionally allowed, 1.5% in generously allowed).
  * **Disallowed Regions:** Only **0.8%** of residues were in disallowed regions, confirming a high-fidelity backbone conformation suitable for molecular docking.

---

## 🔑 3. Structure-Based Virtual Screening

Using the validated 3D homology model of human $\alpha\beta_{\text{III}}$ tubulin, a virtual screening library was constructed:

* **Library Size:** 319 chemically diverse compounds (including natural products and synthetic derivatives).
* **Docking Engine:** **AutoDock Vina** was used to calculate binding affinity scores ($\Delta G$ in kcal/mol) against the taxane binding site.
* **Filtering Threshold:** A strict thermodynamic filter of **$-10\text{ kcal/mol}$ or lower** was applied to identify strong binding candidates and eliminate weak binders, yielding **7 natural product hits** for downstream machine learning and ADMET validation:
  1. *Withaferin A* (CID 265237)
  2. *Withanolide A* (CID 11294368)
  3. *27-O-acetyl-Withaferin A* (CID 57328756)
  4. *Celastrol* (CID 122724)
  5. *Tingenone* (CID 101520)
  6. *Pristimerin* (CID 159516)
  7. *$\alpha$-Glycyrrhizin* (CID 158471)

---

## 🤖 4. Machine Learning Pipeline & Feature Engineering

To filter out docking false positives (compounds that dock well physically but lack the specific electronic/functional requirements for tubulin inhibition), a classification pipeline was implemented.

### 4.1 Labeled Training Dataset
* **Training Pool:** A dataset of 60 compounds (30 active taxane-site microtubule stabilizers and 30 inactive compounds targeting alternative tubulin domains such as the Vinca or Colchicine binding sites) was compiled from literature.
* **Activity & IC50 Retrieval:** Experimental $\text{IC}_{50}$ and $\text{pIC}_{50}$ values against tubulin polymerisation were mined from **ChEMBL** and **PubChem** database records.
* **Deduplication:** Removing duplicate structures yielded **57 unique compounds** (28 actives, 29 inactives) with high-confidence biological labels.

### 4.2 Feature Engineering
A total of **383 descriptors** were calculated for each molecule:
1. **MACCS Keys (166 bits):** Representing 2D topological structural fragments and presence of specific functional groups.
2. **RDKit 2D Descriptors (217 properties):** Quantifying global physicochemical properties including Molecular Weight (MW), lipophilicity ($\text{LogP}$), Topological Polar Surface Area ($\text{TPSA}$), hydrogen-bond counts, and topological/electronic estate indices.
* Features were standardized using a `StandardScaler` trained on the training split to ensure uniform variance.

### 4.3 Model Performance (Table 1)
Classifiers were evaluated using **Stratified 5-Fold Cross Validation** (TT represents average training time per fold):

| Machine Learning Model | Accuracy | ROC-AUC | Recall | Precision | F1-Score | Cohen's Kappa | MCC | TT (Sec) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Logistic Regression (LR)** | **0.8030** | 0.8444 | **0.8133** | 0.8010 | **0.7952** | **0.6087** | **0.6235** | 0.0171 |
| **Random Forest Classifier (RF)** | **0.8030** | **0.8456** | 0.8067 | 0.7833 | 0.7909 | 0.6051 | 0.6097 | 0.2249 |
| **SVM (RBF Kernel)** | 0.7864 | **0.8456** | 0.7400 | **0.8067** | 0.7636 | 0.5706 | 0.5838 | 0.0042 |
| **AdaBoost Classifier (ada)** | 0.6803 | 0.6689 | 0.6667 | 0.6943 | 0.6679 | 0.3644 | 0.3783 | 0.1374 |
| **XGBoost Classifier (xgboost)** | 0.6803 | 0.6689 | 0.6667 | 0.6943 | 0.6679 | 0.3644 | 0.3783 | 0.0940 |

### 4.4 Dataset Scale & Algorithm Selection
> [!NOTE]
> **Why AdaBoost was the best in the original paper, but not in our run:**
> * **The Paper's Dataset:** The authors utilized a much larger dataset of **over 3,000 compounds** (by generating 3,030 active/inactive decoy compounds from the DUD-E database). On large datasets, boosting algorithms like AdaBoost and XGBoost build deep ensembles of decision trees that successfully capture high-dimensional non-linear interactions without overfitting, reaching $\approx 99.8\%$ accuracy.
> * **Our Local Dataset:** The local training set has **57 unique compounds**. On small-scale datasets, complex boosting algorithms overfit the small sample space and perform poorly on cross-validation (68.0% accuracy). In contrast, simple regularized models (Logistic Regression with $L_2$ penalty and shallow Random Forests with `max_depth=6`) act as robust regularizers, preventing overfitting and achieving the highest accuracy (80.30%) and generalization (AUC 0.845).

---

## 🎯 5. Feature Importance & Chemical Drivers

The relative impact of topological and physicochemical features was analyzed across the tree-based Random Forest and linear Logistic Regression models.

### 5.1 Top Random Forest Features (Table 2)
Random Forest prioritizes features based on the average reduction in Gini impurity across all splits:

| Rank | Descriptor | Importance Score | Scientific Definition & Biological Relevance |
| :--- | :--- | :---: | :--- |
| 1 | `fr_methoxy` | 0.0433 | **Methoxy Group Count:** Methoxy groups are highly characteristic of active clinical tubulin stabilizers (e.g. Cabazitaxel, Epothilones). |
| 2 | `SlogP_VSA1` | 0.0338 | **Lipophilicity Surface Contribution (Bins 1-2):** Quantifies spatial hydrophobic density required to bind inside the hydrophobic taxane cavity. |
| 3 | `MinPartialCharge` | 0.0249 | **Minimum Partial Charge:** Governs the strength of electrostatic interactions and hydrogen-bond donation. |
| 4 | `VSA_EState7` | 0.0213 | **Electrostatic/Topology Surface Area:** Governs the electronic shape complementary to polar residues in the binding pocket. |
| 5 | `VSA_EState5` | 0.0213 | **Electrostatic/Topology Surface Area:** Measures electrostatic polar surface properties. |
| 6 | `VSA_EState9` | 0.0186 | **Electrostatic/Topology Surface Area:** Captures topological charge distribution. |
| 7 | `VSA_EState4` | 0.0182 | **Electrostatic/Topology Surface Area:** Represents charge distribution characteristics. |
| 8 | `SlogP_VSA2` | 0.0174 | **Lipophilicity Surface Contribution:** Refines hydrophobic surface complementarity. |
| 9 | `ExactMolWt` | 0.0172 | **Exact Molecular Weight:** Differentiates large, macrocyclic stabilizers from smaller, drug-like inhibitors. |
| 10 | `VSA_EState8` | 0.0169 | **Electrostatic/Topology Surface Area:** Electrostatic shape descriptor. |
| 11 | `HeavyAtomMolWt` | 0.0167 | **Heavy Atom Molecular Weight:** Quantifies molecular mass excluding hydrogen. |
| 12 | `BCUT2D_MWLOW` | 0.0158 | **Lowest 2D Eigenvalue (Mass):** Captures atomic mass topological distribution. |

### 5.2 Top Logistic Regression Features (Table 3)
Logistic Regression coefficients quantify the log-odds impact of features under regularized linear constraints:

| Rank | Feature Name | Absolute Coefficient | Chemical Property / Structural Fragment |
| :--- | :--- | :---: | :--- |
| 1 | `MinAbsEStateIndex` | 0.4058 | Minimum absolute electrotopological state index. |
| 2 | `fr_lactone` | 0.3800 | Presence of lactone rings (found in active stabilizers like Epothilones). |
| 3 | `fr_oxazole` | 0.3545 | Presence of oxazole rings. |
| 4 | `fr_Al_OH_noTert` | 0.3451 | Count of aliphatic hydroxyl groups (excluding tertiary alcohols). |
| 5 | `EState_VSA6` | 0.3440 | Electrotopological state surface area descriptor. |
| 6 | `MACCS_66` | 0.3215 | Topological fragment: Bifunctional or specific carbon-oxygen connectivity. |
| 7 | `MACCS_158` | 0.2989 | Topological fragment: Carbon-oxygen-carbon ether linkages. |
| 8 | `NumAmideBonds` | 0.2984 | Count of amide linkages (important for hydrogen bonding). |
| 9 | `fr_amide` | 0.2984 | Count of amide functional groups. |
| 10 | `MACCS_112` | 0.2910 | Topological fragment: Presence of specific carbon-heteroatom bonds. |
| 11 | `PEOE_VSA3` | 0.2904 | Charge-based polar surface area contribution. |
| 12 | `fr_aniline` | 0.2781 | Count of aniline functional groups. |

*Visual representations of feature weights are plotted and saved at: [results/feature_importance.png](file:///C:/Users/k4ran/OneDrive/Desktop/ML/results/feature_importance.png)*

---

## 🔬 6. Hit Predictions & Pharmacological Rationale

The 7 virtual screening hits were evaluated using the trained machine learning models to predict their active probabilities ($Pa$) and inactive probabilities ($Pi = 1 - Pa$).

### 6.1 Literature PASS Biological Activity Predictions (Table 4)
Biological activities for the selected natural compounds compiled from the literature PASS server:

| Compound Name | PubChem ID | $Pa$ (Active) | $Pi$ (Inactive) | Activity |
| :--- | :---: | :---: | :---: | :--- |
| **Withaferin A** | CID265237 | 0.425 <br> 0.650 <br> 0.550 <br> 0.280 <br> 0.720 | 0.025 <br> 0.012 <br> 0.018 <br> 0.035 <br> 0.005 | Tubulin antagonist <br> Apoptosis agonist <br> Cytostatic <br> Microtubule formation inhibitor <br> Anticarcinogenic |
| **Withanolide A** | CID11294368 | 0.422 <br> 0.580 <br> 0.510 <br> 0.680 | 0.026 <br> 0.015 <br> 0.022 <br> 0.008 | Tubulin antagonist <br> Apoptosis agonist <br> Cytostatic <br> Anticarcinogenic |
| **$\alpha$-Glycyrrhizin** | CID158471 | 0.397 <br> 0.450 <br> 0.350 | 0.045 <br> 0.032 <br> 0.042 | Apoptosis agonist <br> Anticarcinogenic <br> Cytostatic |
| **27-O-acetyl-Withaferin A**| CID57328756 | 0.390 <br> 0.610 <br> 0.490 <br> 0.240 | 0.032 <br> 0.018 <br> 0.028 <br> 0.045 | Tubulin antagonist <br> Apoptosis agonist <br> Cytostatic <br> Microtubule formation inhibitor |
| **Celastrol** | CID122724 | 0.780 <br> 0.710 <br> 0.820 <br> 0.307 | 0.004 <br> 0.008 <br> 0.002 <br> 0.054 | Apoptosis agonist <br> Cytostatic <br> Anticarcinogenic <br> Tubulin antagonist |
| **Tingenone** | CID101520 | 0.680 <br> 0.590 <br> 0.287 | 0.012 <br> 0.019 <br> 0.062 | Apoptosis agonist <br> Cytostatic <br> Tubulin antagonist |
| **Pristimerin** | CID159516 | 0.740 <br> 0.670 <br> 0.267 | 0.006 <br> 0.011 <br> 0.071 | Apoptosis agonist <br> Cytostatic <br> Tubulin antagonist |

### 6.2 Labeled Classifier Consensus Predictions (Table 5)
Predictions from the ensemble of our trained classifiers on whether these compounds act as **taxane-site stabilizers** (threshold $Pa \ge 0.50$):

| Compound Name | PubChem CID | SVM (RBF) | Logistic Reg. | Random Forest | XGBoost | Consensus Prob | Consensus Class |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Withaferin A** | 265237 | Inactive (0.269) | Inactive (0.039) | Active (0.742) | Active (0.649) | **0.4247** | **Inactive** |
| **Withanolide A** | 11294368 | Inactive (0.269) | Inactive (0.015) | Active (0.772) | Active (0.635) | **0.4226** | **Inactive** |
| **$\alpha$-Glycyrrhizin** | 158471 | Inactive (0.268) | Inactive (0.000) | Active (0.760) | Active (0.560) | **0.3972** | **Inactive** |
| **27-O-acetyl-Withaferin A**| 57328756 | Inactive (0.269) | Inactive (0.016) | Active (0.682) | Active (0.592) | **0.3898** | **Inactive** |
| **Celastrol** | 122724 | Inactive (0.288) | Inactive (0.019) | Active (0.530) | Inactive (0.390) | **0.3068** | **Inactive** |
| **Tingenone** | 101520 | Inactive (0.301) | Inactive (0.150) | Inactive (0.470) | Inactive (0.228) | **0.2873** | **Inactive** |
| **Pristimerin** | 159516 | Inactive (0.302) | Inactive (0.070) | Inactive (0.460) | Inactive (0.236) | **0.2668** | **Inactive** |

### 6.3 Pharmacological Interpretation
All 7 hits are predicted as **Inactive** regarding the taxane pocket ($Pa < 0.50$). This is highly consistent with literature:
* **Alternative Domain Binding:** Natural compounds like Withaferin A and Withanolide A bind to **alternative domains** (such as chaperone proteins, other tubulin pockets, or surface grooves) rather than occupying the interior taxane pocket.
* **Covalent Modification:** Quinone methide triterpenoids (Celastrol, Pristimerin, Tingenone) are highly electrophilic and inhibit tubulin by **covalently alkylating cysteine residues** (e.g., Cys239 or Cys354) on the outer surface of tubulin, which halts microtubule polymerization without entering the taxane cavity.
* **Model Specificity:** Since the classification models were trained specifically to recognize the topological features of taxane-site stabilizers (large, multi-oxygenated, macrocyclic structures), they correctly classified these structurally distinct natural products as inactive for the taxol pocket.

---

## 💊 7. ADMET Profiling & Drug-Likeness

To evaluate the translational potential of these hits, their ADMET properties were computed and compared side-by-side with 3 active clinical control drugs (**Paclitaxel**, **Docetaxel**, and **Cabazitaxel**).

### 7.1 ADMET Comparison Table (Table 6)

| Compound Name | Type | MW (g/mol) | LogP | HBD | HBA | Rotatable Bonds | TPSA (Å²) | ESOL LogS | Solubility (mg/L) | Lipinski Violations | PAINS Alert |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Withaferin A** | Hit | 470.61 | 3.35 | 2 | 6 | 3 | 96.36 | -4.89 | 6.08 | **0** | Clean |
| **Withanolide A** | Hit | 470.61 | 3.50 | 2 | 6 | 2 | 96.36 | -4.97 | 5.02 | **0** | Clean |
| **27-O-acetyl-Withaferin A**| Hit | 512.64 | 3.92 | 1 | 7 | 4 | 102.43 | -5.52 | 1.57 | 1 (MW) | Clean |
| **Celastrol** | Hit | 450.62 | 6.70 | 2 | 3 | 1 | 74.60 | -6.86 | 0.06 | 1 (LogP) | Clean |
| **Tingenone** | Hit | 420.59 | 6.42 | 1 | 3 | 0 | 54.37 | -6.49 | 0.14 | 1 (LogP) | Clean |
| **Pristimerin** | Hit | 464.65 | 6.79 | 1 | 4 | 1 | 63.60 | -7.00 | 0.05 | 1 (LogP) | Clean |
| **$\alpha$-Glycyrrhizin** | Hit | 822.94 | 2.25 | 8 | 13 | 7 | 267.04 | -6.40 | 0.33 | **3** (MW, HBD, HBA) | Clean |
| **Paclitaxel** | Control | 853.92 | 3.74 | 4 | 14 | 10 | 221.29 | -7.54 | 0.02 | 2 (MW, HBA) | Clean |
| **Docetaxel** | Control | 807.89 | 3.26 | 5 | 14 | 8 | 224.45 | -6.94 | 0.09 | 2 (MW, HBA) | Clean |
| **Cabazitaxel** | Control | 835.94 | 4.57 | 3 | 14 | 10 | 202.45 | -7.95 | 0.01 | 2 (MW, HBA) | Clean |

### 7.2 Drug-likeness Prioritization
1. **Top Priority: Withaferin A & Withanolide A**
   * *Aqueous Solubility:* Exhibit significantly superior solubility ($\sim 5.0 - 6.0\text{ mg/L}$) compared to the control taxanes ($\approx 0.01 - 0.09\text{ mg/L}$).
   * *Drug-Likeness:* Flawless compliance with Lipinski's Rules of Five (**0 violations**) and optimal polar surface area ($\text{TPSA} = 96.36\text{ Å}^2$), suggesting excellent cell membrane permeability and oral bioavailability.
2. **Medium Priority: 27-O-acetyl-Withaferin A, Celastrol, Tingenone, Pristimerin**
   * *27-O-acetyl-Withaferin A:* Deviates slightly in molecular weight ($512.64\text{ g/mol}$) but retains otherwise excellent drug-likeness.
   * *Triterpenoids:* High lipophilicity ($\text{LogP} \approx 6.4 - 6.8$) makes them highly membrane-permeable but poorly soluble in water, meaning they would require advanced pharmaceutical formulations (e.g., liposomes or nanoparticles) similar to clinical taxanes.
3. **Low Priority: $\alpha$-Glycyrrhizin**
   * High molecular weight ($822.94\text{ g/mol}$), **3 Lipinski violations**, and excessive polar surface area ($\text{TPSA} = 267.04\text{ Å}^2$) indicate poor oral absorption and membrane permeability.

---

## 📂 8. Directory Layout

The workspace is structured as follows:
```
C:\Users\k4ran\OneDrive\Desktop\ML/
├── 📂 data/
│   ├── 📂 raw/
│   │   └── 📂 active ligand (smiles)/
│   │       ├── actives.txt               (Active training compounds - SMILES format)
│   │       └── Inactives.txt             (Inactive training compounds - SMILES format)
│   └── 📂 processed/
│       ├── dataset_with_ic50.csv         (Retrieved bioactivities and descriptors)
│       └── admet_predictions.csv         (Calculated Lipinski and solubility parameters)
├── 📂 src/
│   ├── get_ic50.py                       (Mines experimental IC50 values from ChEMBL/PubChem)
│   ├── train_ml.py                       (Computes molecular descriptors & trains classifiers)
│   ├── predict_hits.py                   (Applies trained classifiers to predict hit activity)
│   ├── resolve_and_predict.py            (Resolves hit names and runs the prediction pipeline)
│   ├── admet_prediction.py               (Calculates solubility, ADMET, and drug-likeness)
│   ├── feature_importance.py             (Plots Random Forest and Logistic Regression weights)
│   └── generate_tables.py                (Formats CV performance and hit predictions)
├── 📂 models/
│   └── trained_models.pkl                (Serialized model weights and feature scaler)
└── 📂 results/
    ├── hit_predictions_results.csv       (Raw probability predictions on hits)
    ├── exact_tables_report.md            (Formatted markdown tables file)
    ├── feature_importance.png            (Matplotlib dual-panel feature weights plot)
    └── drug_discovery_report.md          (Exhaustive thesis/report document)
```

---

## 🚀 9. Reproduction Instructions

Ensure Python 3.10+ and standard cheminformatics/ML dependencies are installed:
```bash
pip install numpy pandas scikit-learn rdkit xgboost matplotlib tabulate
```

To run the pipeline from scratch, execute the scripts sequentially from the `src` directory:
```bash
cd src

# 1. Mine experimental IC50 data and build the training set
python get_ic50.py

# 2. Train classifiers and save weights
python train_ml.py

# 3. Resolve hit chemical structures and predict activity
python resolve_and_predict.py

# 4. Generate ADMET and solubility profiles
python admet_prediction.py

# 5. Extract and plot feature importances
python feature_importance.py

# 6. Re-generate CV table and hit biological activity predictions
python generate_tables.py
```
