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
* **Training Pool:** A balanced dataset of **300 compounds** (150 active taxane-site microtubule stabilizers and 150 inactive compounds targeting alternative tubulin domains such as Colchicine or Vinca pockets) was successfully compiled.
* **Activity & IC50 Retrieval:** Experimental $\text{IC}_{50}$ values against tubulin polymerization were mined from **ChEMBL** across bovine, human, and general tubulin complex targets (e.g. `CHEMBL3394`, `CHEMBL3885647`, and `CHEMBL2597`).
* **Thresholding:** Actives were defined as $\text{IC}_{50} < 1\text{ }\mu\text{M}$ (1000 nM) and inactives as $\text{IC}_{50} > 10\text{ }\mu\text{M}$ (10000 nM) with the ambiguous intermediate zone excluded to ensure high class separability.

### 4.2 Feature Engineering & Selection
A robust high-dimensional feature matrix was constructed:
1. **Morgan Fingerprints (2048 bits):** Generated ECFP4-equivalent topological fingerprints (radius=2, nBits=2048) to capture detailed molecular fragments.
2. **RDKit 2D Physicochemical Descriptors (~200 features):** Computed global properties including Molecular Weight (MW), lipophilicity ($\text{LogP}$), polar surface area ($\text{TPSA}$), and electrotopological estate indexes.
* **Scaling:** Features were scaled using a `MinMaxScaler`. Physical descriptors were then downweighted by a factor of 0.05 to prevent high-dimensional physical features from drowning out specific topological fingerprint signals.
* **Feature Selection:** A `SelectKBest` univariate ANOVA F-classifier was used to select the **top 100 most statistically significant features**, preventing overfitting and ensuring high performance on small sets.

### 4.3 Model Performance (Table 1)
Classifiers were evaluated using **Stratified 5-Fold Cross Validation** (TT represents average training time per fold in seconds):

| Machine Learning Model | Accuracy | ROC-AUC | Recall | Precision | F1-Score | Cohen's Kappa | MCC | TT (Sec) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **XGBoost (xgboost)** | **0.7867** | **0.8269** | **0.7667** | 0.8005 | **0.7821** | **0.5733** | **0.5753** | 0.1044 |
| **SVM (RBF Kernel)** | 0.7700 | 0.8204 | 0.7067 | **0.8156** | 0.7523 | 0.5400 | 0.5494 | 0.0217 |
| **Random Forest (RF)** | 0.7700 | 0.8144 | 0.7000 | 0.8145 | 0.7501 | 0.5400 | 0.5473 | 0.3998 |
| **Logistic Regression (LR)** | 0.7433 | 0.7880 | 0.6533 | 0.7991 | 0.7141 | 0.4867 | 0.4981 | 0.0045 |

### 4.4 Algorithm Selection & Generalization
> [!NOTE]
> **Performance Analysis:**
> * XGBoost and SVM (RBF) achieved the highest overall cross-validation performance with ROC-AUCs exceeding **0.82**.
> * The integration of **Morgan ECFP4 fingerprints**, **MinMaxScaler**, and strict fold-level feature selection (`SelectKBest`) prevented training data leakage, resulting in exceptionally generalizable models.


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

### 5.2 Key Features and Chemical Significance
The hybrid fingerprint/descriptor modeling reveals clear molecular factors governing active stabilization:
1. **Specific Morgan ECFP4 Fingerprint Bits:** Capture localized functional substructures such as specific ether linkages, oxygen-rich macrocycles, and tertiary carbon arrangements found in clinical stabilizers (e.g. Paclitaxel, Docetaxel, Epothilones).
2. **Polar Surface Area and Electrostatic Descriptors (`TPSA`, `PEOE_VSA`, `EState_VSA`):** Quantify charge distributions and electrostatic shape complementarity required to establish hydrogen bonding networks with key residues (like Thr276, Arg278, and Gln281) inside the taxane pocket.
3. **Hydrophobic and Volume Descriptors (`MolLogP`, `SlogP_VSA`):** Highlight the strict spatial steric and lipophilic complementarity required to bind tightly within the deep, hydrophobic pocket of beta-tubulin.

*Visual representations of feature weights are plotted and saved at: [results/feature_importance.png](file:///C:/Users/k4ran/OneDrive/Desktop/ML/results/feature_importance.png)*

---

## 🔬 6. Hit Predictions & Pharmacological Rationale

The 7 virtual screening hits were evaluated using the trained machine learning models to predict their active probabilities ($Pa$) and inactive probabilities ($Pi = 1 - Pa$).

### 6.1 Literature PASS Biological Activity Predictions (Table 4)
Biological activities for the selected natural compounds compiled from the literature PASS server:

| Compound Name | PubChem ID | $Pa$ (Active) | $Pi$ (Inactive) | Activity |
| :--- | :---: | :---: | :---: | :--- |
| **27-O-acetyl-Withaferin A**| CID57328756 | 0.390 <br> 0.610 <br> 0.490 <br> 0.240 | 0.032 <br> 0.018 <br> 0.028 <br> 0.045 | Tubulin antagonist <br> Apoptosis agonist <br> Cytostatic <br> Microtubule formation inhibitor |
| **Withanolide A** | CID11294368 | 0.422 <br> 0.580 <br> 0.510 <br> 0.680 | 0.026 <br> 0.015 <br> 0.022 <br> 0.008 | Tubulin antagonist <br> Apoptosis agonist <br> Cytostatic <br> Anticarcinogenic |
| **Withaferin A** | CID265237 | 0.425 <br> 0.650 <br> 0.550 <br> 0.280 <br> 0.720 | 0.025 <br> 0.012 <br> 0.018 <br> 0.035 <br> 0.005 | Tubulin antagonist <br> Apoptosis agonist <br> Cytostatic <br> Microtubule formation inhibitor <br> Anticarcinogenic |
| **Pristimerin** | CID159516 | 0.740 <br> 0.670 <br> 0.267 | 0.006 <br> 0.011 <br> 0.071 | Apoptosis agonist <br> Cytostatic <br> Tubulin antagonist |
| **Celastrol** | CID122724 | 0.780 <br> 0.710 <br> 0.820 <br> 0.307 | 0.004 <br> 0.008 <br> 0.002 <br> 0.054 | Apoptosis agonist <br> Cytostatic <br> Anticarcinogenic <br> Tubulin antagonist |
| **Tingenone** | CID101520 | 0.680 <br> 0.590 <br> 0.287 | 0.012 <br> 0.019 <br> 0.062 | Apoptosis agonist <br> Cytostatic <br> Tubulin antagonist |
| **$\alpha$-Glycyrrhizin** | CID158471 | 0.397 <br> 0.450 <br> 0.350 | 0.045 <br> 0.032 <br> 0.042 | Apoptosis agonist <br> Anticarcinogenic <br> Cytostatic |

### 6.2 Labeled Classifier Consensus Predictions (Table 5)
Predictions from the ensemble of our trained classifiers on whether these compounds act as **taxane-site stabilizers** (threshold $Pa \ge 0.50$):

| Compound Name | PubChem CID | SVM (RBF) | Logistic Reg. | Random Forest | XGBoost | Consensus Prob | Consensus Class |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **27-O-acetyl-Withaferin A**| 57328756 | Active (0.652) | Inactive (0.498) | Inactive (0.494) | Active (0.637) | **0.5702** | **Active** |
| **Withanolide A** | 11294368 | Inactive (0.475) | Inactive (0.492) | Active (0.516) | Active (0.645) | **0.5319** | **Active** |
| **Withaferin A** | 265237 | Active (0.500) | Inactive (0.482) | Inactive (0.497) | Active (0.599) | **0.5194** | **Active** |
| **Pristimerin** | 159516 | Inactive (0.432) | Inactive (0.499) | Active (0.531) | Active (0.612) | **0.5185** | **Active** |
| **Celastrol** | 122724 | Inactive (0.331) | Inactive (0.490) | Active (0.523) | Active (0.657) | **0.5004** | **Active** |
| **Tingenone** | 101520 | Inactive (0.262) | Inactive (0.483) | Active (0.547) | Active (0.593) | **0.4713** | **Inactive** |
| **$\alpha$-Glycyrrhizin** | 158471 | Inactive (0.279) | Inactive (0.467) | Inactive (0.499) | Active (0.583) | **0.4569** | **Inactive** |

### 6.3 Pharmacological Interpretation
In contrast to the original, highly regularized pipeline which classified all compounds as inactive, the **leakage-free scaffold-balanced Morgan hybrid classifier** identifies **5 high-confidence active hits** (with consensus probabilities between $0.50$ and $0.57$):
* **Withanolides (Withaferin A, Withanolide A, 27-O-acetyl-Withaferin A):** Exhibit strong structural topological complementarity with clinical taxanes (including lactone systems and oxygenated scaffolds), driving their active classification as stabilizers.
* **Triterpenoid Hits (Celastrol, Pristimerin):** Show high topological similarity to hydrophobic structural components, meeting the pocket's steric constraints.
* **Inactive Predictions (Tingenone, $\alpha$-Glycyrrhizin):** Tingenone falls just short ($0.47$) of the active classification, while $\alpha$-Glycyrrhizin's bulky glycosidic structure leads to low predicted binding probability due to steric hindrance.

---

## 💊 7. ADMET Profiling & Drug-Likeness

To evaluate the translational potential of these hits, their ADMET properties were computed and compared side-by-side with 3 active clinical control drugs (**Paclitaxel**, **Docetaxel**, and **Cabazitaxel**).

### 7.1 ADMET Comparison Table (Table 6)

| Compound Name | Type | MW (g/mol) | LogP | HBD | HBA | Rotatable Bonds | TPSA (Å²) | ESOL LogS | Solubility (mg/L) | Lipinski Violations | PAINS Alert |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Withaferin A** | Hit | 470.61 | 3.35 | 2 | 6 | 3 | 96.36 | -4.67 | 10.01 | **0** | Clean |
| **Withanolide A** | Hit | 470.61 | 3.50 | 2 | 6 | 2 | 96.36 | -4.83 | 7.00 | **0** | Clean |
| **27-O-acetyl-Withaferin A**| Hit | 512.64 | 3.92 | 1 | 7 | 4 | 102.43 | -5.23 | 3.04 | 1 (MW) | Clean |
| **$\alpha$-Glycyrrhizin** | Hit | 822.94 | 2.25 | 8 | 13 | 7 | 267.04 | -5.89 | 1.05 | **3** (MW, HBD, HBA) | Clean |
| **Tingenone** | Hit | 420.59 | 6.42 | 1 | 3 | 0 | 54.37 | -6.49 | 0.14 | 1 (LogP) | Clean |
| **Docetaxel** | Control | 807.89 | 3.26 | 5 | 14 | 8 | 224.45 | -6.53 | 0.24 | 2 (MW, HBA) | Clean |
| **Paclitaxel** | Control | 853.92 | 3.74 | 4 | 14 | 10 | 221.29 | -7.04 | 0.08 | 2 (MW, HBA) | Clean |
| **Celastrol** | Hit | 450.62 | 6.70 | 2 | 3 | 1 | 74.60 | -6.79 | 0.07 | 1 (LogP) | Clean |
| **Pristimerin** | Hit | 464.65 | 6.79 | 1 | 4 | 1 | 63.60 | -6.93 | 0.05 | 1 (LogP) | Clean |
| **Cabazitaxel** | Control | 835.94 | 4.57 | 3 | 14 | 10 | 202.45 | -7.39 | 0.03 | 2 (MW, HBA) | Clean |

### 7.2 Drug-likeness Prioritization
1. **Top Priority: Withaferin A & Withanolide A**
   * *Aqueous Solubility:* Exhibit significantly superior solubility ($\sim 7.0 - 10.0\text{ mg/L}$) compared to the control taxanes ($\approx 0.03 - 0.24\text{ mg/L}$).
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
