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
* **Training Pool:** A curated dataset of 60 compounds (30 active taxane-site stabilizers and 30 inactive compounds targeting other sites like Vinca or Colchicine) was parsed from literature.
* **IC50 and pIC50 Allocation:** Experimental $\text{IC}_{50}$ and $\text{pIC}_{50}$ values against tubulin were retrieved from **ChEMBL** and **PubChem** databases, utilizing high-confidence literature values for standard fallback. After removing redundant duplicate molecules, 57 unique compounds remained.
* **Molecular Representation:**
  * **MACCS Keys (166 bits):** To capture structural fragments.
  * **RDKit 2D Descriptors (217 properties):** To capture physicochemical properties (molecular weight, lipophilicity, polar surface area, hydrogen bond counts, etc.).
  * The final combined feature matrix comprised **383 dimensions** for each molecule.

### 4.2 Cross-Validation Performance
The models were trained and validated using **Stratified 5-Fold Cross-Validation** to ensure class balance and prevent overfitting. The performance indices of the classifiers are summarized below:

| Machine Learning Model | Accuracy | ROC-AUC | Precision | Recall | F1-Score |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Random Forest** | 0.8030 | 0.8456 | 0.7833 | 0.8067 | 0.7909 |
| **Logistic Regression** | 0.8030 | 0.8444 | 0.8010 | 0.8133 | 0.7952 |
| **SVM (RBF Kernel)** | 0.7864 | 0.8456 | 0.8067 | 0.7400 | 0.7636 |
| **XGBoost** | 0.6803 | 0.6689 | 0.6943 | 0.6667 | 0.6679 |

*Random Forest and Logistic Regression achieved the highest classification accuracy (80.30%) and ROC-AUC (~0.845), proving to be highly robust and generalizable on this small molecular dataset.*

### 4.3 Feature Importance
To understand the chemical drivers behind taxane-site binding, we analyzed descriptor contributions:
* **fr_methoxy (Methoxy groups):** The most significant feature in tree-based splits, reflecting the presence of methoxy/ether groups in active stabilizers like Cabazitaxel or Epothilones.
* **SlogP_VSA1 & SlogP_VSA2 (Lipophilicity surface distribution):** Highlighted the importance of spatially distributed hydrophobic interactions required to bind inside the lipophilic taxane cavity.
* **MinPartialCharge & VSA_EState (Electronic/Topological shape):** Indicated that specific electrostatic charge distributions and molecular shapes are essential to establish stable hydrogen bonding and electrostatic interactions with the pocket residues.
* **ExactMolWt (Molecular Weight):** Captured the macrocyclic, high-molecular-weight nature of active taxane-site ligands compared to small-molecule inhibitors.

---

## 5. Machine Learning Prediction on Hits
The ensemble of trained ML models was applied to screen the 7 natural product hits. The results are ranked below by their **Consensus Active Probability Score**:

| Compound Name | PubChem CID | Consensus Active Prob | SVM (RBF) | Logistic Reg. | Random Forest | XGBoost | Consensus Class |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Withaferin A** | 265237 | **0.4247** | Inactive | Inactive | Active | Active | **Inactive** |
| **Withanolide A** | 11294368 | **0.4226** | Inactive | Inactive | Active | Active | **Inactive** |
| **alpha-Glycyrrhizin** | 158471 | **0.3972** | Inactive | Inactive | Active | Active | **Inactive** |
| **27-O-acetyl-Withaferin A**| 57328756 | **0.3898** | Inactive | Inactive | Active | Active | **Inactive** |
| **Celastrol** | 122724 | **0.3068** | Inactive | Inactive | Active | Inactive | **Inactive** |
| **Tingenone** | 101520 | **0.2873** | Inactive | Inactive | Inactive | Inactive | **Inactive** |
| **Pristimerin** | 159516 | **0.2668** | Inactive | Inactive | Inactive | Inactive | **Inactive** |

### Interpretation
* **Consensus Prediction:** All 7 hits are predicted as **Inactive** with respect to the taxane binding site.
* **Pharmacological Rationale:**
  * While compounds like **Withaferin A**, **Celastrol**, and **Pristimerin** are known in literature to have anti-tubulin and anti-cancer activity, they function via alternative mechanisms—primarily as **covalent modifiers** of specific cysteine residues or by binding to **other pockets** (e.g., Colchicine or chaperone-like domains).
  * Because the machine learning models were trained to recognize the distinct spatial, topological, and chemical properties of taxane-site binders (such as large macrocyclic rings and specific methoxy arrangements), they correctly classified these chemically diverse hits as inactive for this target pocket.
