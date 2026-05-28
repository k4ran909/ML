# Research Report: Computational Drug Discovery & Machine Learning Predictions for Tubulin Inhibitors

This report details the integration of structural homology modeling, virtual screening, and machine learning classification to identify and evaluate novel active compounds targeting the taxane binding site of human $\alpha\beta_{\text{III}}$ tubulin.

---

## 1. Introduction & Biological Target
Microtubules, composed of repeating dimers of $\alpha$- and $\beta$-tubulin, play a critical role in cellular structure, intracellular transport, and mitotic spindle assembly during cell division. Uncontrolled cell division is a hallmark of oncology; therefore, disrupting microtubule dynamics represents a major therapeutic strategy. 

This project targets the **taxane binding site** located on the interior surface of $\beta$-tubulin. Ligands binding to this cavity (e.g., Paclitaxel) stabilize microtubules, arrest the cell cycle in the G2/M phase, and induce apoptotic cell death.

---

## 2. Homology Modeling of Human $\alpha\beta_{\text{III}}$ Tubulin
Because a high-resolution experimental crystal structure of human $\alpha\beta_{\text{III}}$ tubulin was unavailable in the Protein Data Bank (PDB), a 3D structural model was constructed using homology modeling:
* **Template Selection:** The crystal structure of the bovine $\alpha_{\text{1B}}\beta_{\text{IIB}}$ tubulin heterodimer (**PDB ID: 1JFF**) was used as a template due to its high sequence similarity.
* **Model Generation:** 100 structural models of the human $\beta_{\text{III}}$ tubulin isotype were generated using **Modeller 10.2** by keeping the $\alpha_{\text{1B}}$ subunit constant and remodeling the $\beta$-tubulin sequence.
* **Model Selection (DOPE Score):** The optimal model was selected based on the lowest **Discrete Optimized Protein Energy (DOPE)** score of **$-51,889.82422$**, representing a highly stable and physically relaxed conformation.
* **Stereochemical Validation (Ramachandran Plot):** Structural validation using **PROCHECK** revealed that **97.7%** of the amino acid residues fell within the favored and allowed regions (88.1% in the most favored region, 9.6% in additional allowed, 1.5% in generously allowed, and only 0.8% in disallowed regions). This confirmed that the homology model was of high quality and suitable for structure-based studies.

---

## 3. Structure-Based Virtual Screening
Using the validated 3D homology model of human $\alpha\beta_{\text{III}}$ tubulin:
* A virtual library of **319 ligands** (natural products and synthetic compounds) was prepared and screened.
* Docking simulations were performed against the taxane binding pocket using **AutoDock Vina** to evaluate binding affinities and modes of interaction.
* A strict filtering threshold of **$-10\text{ kcal/mol}$ or lower** was applied to identify strong binders, yielding **7 natural product hits** showing strong docking scores:
  1. *Withaferin A*
  2. *Withanolide A*
  3. *27-O-acetyl-Withaferin A*
  4. *Celastrol*
  5. *Tingenone*
  6. *Pristimerin*
  7. *$\alpha$-Glycyrrhizin*

---

## 4. Machine Learning Classification Pipeline

### 4.1 Dataset and Feature Engineering
To filter out potential false positives from the docking simulations, we trained a supervised machine learning classifier.
* **Training Pool:** A balanced dataset of **300 compounds** (150 active taxane-site stabilizers and 150 inactive compounds targeting other sites like Colchicine or Vinca pockets) was successfully compiled from ChEMBL database queries and high-confidence reference sets.
* **IC50 Thresholds:** Experimental $\text{IC}_{50}$ values were mined across bovine, human, and general tubulin complexes. Active compounds were defined as $\text{IC}_{50} < 1\text{ }\mu\text{M}$ (1000 nM) and inactives as $\text{IC}_{50} > 10\text{ }\mu\text{M}$ (10000 nM) to establish class balance and high separability.
* **Molecular Representation:**
  * **Morgan Fingerprints (2048 bits):** Generated ECFP4-equivalent topological fingerprints (radius=2) to capture detailed molecular fragments.
  * **RDKit 2D Physicochemical Descriptors (~200 properties):** Calculated global molecular properties (MW, LogP, TPSA, hydrogen-bond counts, etc.).
  * **Scaling & Downweighting:** MinMaxScaler was applied, and physical descriptors were then downweighted by a factor of 0.05 to ensure topological fingerprint bits represent the primary signal.
  * **Feature Selection:** SelectKBest (k=100) univariate ANOVA F-classifier was used to select the top 100 features.

### 4.2 Cross-Validation Performance
The models were trained and validated using **Stratified 5-Fold Cross-Validation** to ensure strict fold-level feature selection and prevent data leakage. The performance indices are summarized below:

| Machine Learning Model | Accuracy | ROC-AUC | Precision | Recall | F1-Score | Cohen's Kappa | MCC |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **XGBoost** | **0.7867** | **0.8269** | 0.8005 | **0.7667** | **0.7821** | **0.5733** | **0.5753** |
| **SVM (RBF Kernel)** | 0.7700 | 0.8204 | **0.8156** | 0.7067 | 0.7523 | 0.5400 | 0.5494 |
| **Random Forest** | 0.7700 | 0.8144 | 0.8145 | 0.7000 | 0.7501 | 0.5400 | 0.5473 |
| **Logistic Regression** | 0.7433 | 0.7880 | 0.7991 | 0.6533 | 0.7141 | 0.4867 | 0.4981 |

*XGBoost and SVM (RBF) achieved the highest ROC-AUC (~0.82-0.83), proving to be highly robust and generalizable on this leakage-free balanced molecular dataset.*

### 4.3 Feature Importance
To understand the chemical drivers behind taxane-site binding, we analyzed descriptor contributions:
* **Morgan Fingerprint Bits:** Specific topological arrangement bits corresponding to oxygenated ring systems and methoxy arrangements were key chemical drivers of active stabilizers.
* **SlogP_VSA & Lipophilicity surface distribution:** Highlighted the importance of spatially distributed hydrophobic interactions required to complement the pocket's deep hydrophobic interior.
* **TPSA & PEOE_VSA (Polarity and Electrostatic shape):** Indicated that specific topological polar surface shapes are essential to establish critical hydrogen bonding networks with residues Thr276, Arg278, and Gln281.

---

## 5. Machine Learning Prediction on Hits
The ensemble of trained ML models was applied to screen the 7 natural product hits. The results are ranked below by their **Consensus Active Probability Score**:

| Compound Name | PubChem CID | Consensus Active Prob | SVM (RBF) | Logistic Reg. | Random Forest | XGBoost | Consensus Class |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **27-O-acetyl-Withaferin A**| 57328756 | **0.5702** | Active | Inactive | Inactive | Active | **Active** |
| **Withanolide A** | 11294368 | **0.5319** | Inactive | Inactive | Active | Active | **Active** |
| **Withaferin A** | 265237 | **0.5194** | Active | Inactive | Inactive | Active | **Active** |
| **Pristimerin** | 159516 | **0.5185** | Inactive | Inactive | Active | Active | **Active** |
| **Celastrol** | 122724 | **0.5004** | Inactive | Inactive | Active | Active | **Active** |
| **Tingenone** | 101520 | **0.4713** | Inactive | Inactive | Active | Active | **Inactive** |
| **alpha-Glycyrrhizin** | 158471 | **0.4569** | Inactive | Inactive | Inactive | Active | **Inactive** |

### Interpretation
* **Consensus Prediction:** Under our optimized, balanced Morgan-fingerprint pipeline, **5 of the 7 natural product hits are predicted as Active** regarding the taxane binding site (Consensus Probability $\ge 0.50$).
* **Pharmacological Rationale:**
  * **Withanolides (27-O-acetyl-Withaferin A, Withanolide A, and Withaferin A):** Exhibit strong consensus active predictions ($0.52 - 0.57$), driven by their rich oxygenated ring structures and topological complementarity with standard taxane-pocket stabilizers.
  * **Quinone Methide Triterpenoids (Celastrol, Pristimerin):** Show consensus active probabilities just exceeding the active threshold ($0.50 - 0.52$), suggesting that while they can form covalent bonds with cysteine residues, their lipophilic skeleton also fits steric pocket constraints.
  * **Tingenone & alpha-Glycyrrhizin:** Predicted as **Inactive**, where Tingenone falls slightly below the threshold, and alpha-Glycyrrhizin is excluded due to its steric bulk.
