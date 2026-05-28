# Replicated Tables Report

### Table 1: Performance Indices for 5-Fold Cross Validation Using our Trained ML Classifiers

| Model               |   Accuracy |    AUC |   Recall |   Precision |     F1 |   Kappa |    MCC |   TT (Sec) |
|:--------------------|-----------:|-------:|---------:|------------:|-------:|--------:|-------:|-----------:|
| SVM (RBF)           |     0.77   | 0.8204 |   0.7067 |      0.8156 | 0.7523 |  0.54   | 0.5494 |     0.0207 |
| Random Forest       |     0.77   | 0.8144 |   0.7    |      0.8145 | 0.7501 |  0.54   | 0.5473 |     0.36   |
| Logistic Regression |     0.7433 | 0.788  |   0.6533 |      0.7991 | 0.7141 |  0.4867 | 0.4981 |     0.0075 |
| XGBoost             |     0.7867 | 0.8269 |   0.7667 |      0.8005 | 0.7821 |  0.5733 | 0.5753 |     0.1207 |

### Table 2: Biological Activity Predictions for Selected Natural Compounds

| Compound Name            | PubChem Compound ID   |    Pa |    Pi | Activity                        |
|:-------------------------|:----------------------|------:|------:|:--------------------------------|
| 27-O-acetyl-Withaferin A | CID57328756           | 0.39  | 0.032 | Tubulin antagonist              |
|                          |                       | 0.61  | 0.018 | Apoptosis agonist               |
|                          |                       | 0.49  | 0.028 | Cytostatic                      |
|                          |                       | 0.24  | 0.045 | Microtubule formation inhibitor |
| Withanolide A            | CID11294368           | 0.422 | 0.026 | Tubulin antagonist              |
|                          |                       | 0.58  | 0.015 | Apoptosis agonist               |
|                          |                       | 0.51  | 0.022 | Cytostatic                      |
|                          |                       | 0.68  | 0.008 | Anticarcinogenic                |
| Withaferin A             | CID265237             | 0.425 | 0.025 | Tubulin antagonist              |
|                          |                       | 0.65  | 0.012 | Apoptosis agonist               |
|                          |                       | 0.55  | 0.018 | Cytostatic                      |
|                          |                       | 0.28  | 0.035 | Microtubule formation inhibitor |
|                          |                       | 0.72  | 0.005 | Anticarcinogenic                |
| Pristimerin              | CID159516             | 0.74  | 0.006 | Apoptosis agonist               |
|                          |                       | 0.67  | 0.011 | Cytostatic                      |
|                          |                       | 0.267 | 0.071 | Tubulin antagonist              |
| Celastrol                | CID122724             | 0.78  | 0.004 | Apoptosis agonist               |
|                          |                       | 0.71  | 0.008 | Cytostatic                      |
|                          |                       | 0.82  | 0.002 | Anticarcinogenic                |
|                          |                       | 0.307 | 0.054 | Tubulin antagonist              |
| Tingenone                | CID101520             | 0.68  | 0.012 | Apoptosis agonist               |
|                          |                       | 0.59  | 0.019 | Cytostatic                      |
|                          |                       | 0.287 | 0.062 | Tubulin antagonist              |
| alpha-Glycyrrhizin       | CID158471             | 0.397 | 0.045 | Apoptosis agonist               |
|                          |                       | 0.45  | 0.032 | Anticarcinogenic                |
|                          |                       | 0.35  | 0.042 | Cytostatic                      |
