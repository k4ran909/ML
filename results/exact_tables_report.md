# Replicated Tables Report

### Table 1: Performance Indices for 5-Fold Cross Validation Using our Trained ML Classifiers

| Model                               |   Accuracy |    AUC |   Recall |   Precision |     F1 |   Kappa |    MCC |   TT (Sec) |
|:------------------------------------|-----------:|-------:|---------:|------------:|-------:|--------:|-------:|-----------:|
| AdaBoost Classifier (ada)           |     0.6803 | 0.6689 |   0.6667 |      0.6943 | 0.6679 |  0.3644 | 0.3783 |     0.0817 |
| Extreme Gradient Boosting (xgboost) |     0.6803 | 0.6689 |   0.6667 |      0.6943 | 0.6679 |  0.3644 | 0.3783 |     0.0595 |
| Logistic Regression (LR)            |     0.803  | 0.8444 |   0.8133 |      0.801  | 0.7952 |  0.6087 | 0.6235 |     0.0104 |
| Random Forest Classifier (RF)       |     0.803  | 0.8456 |   0.8067 |      0.7833 | 0.7909 |  0.6051 | 0.6097 |     0.1153 |
| SVM - RBF Kernel                    |     0.7864 | 0.8456 |   0.74   |      0.8067 | 0.7636 |  0.5706 | 0.5838 |     0.0023 |

### Table 2: Biological Activity Predictions for Selected Natural Compounds

| Compound Name            | PubChem Compound ID   |    Pa |    Pi | Activity                        |
|:-------------------------|:----------------------|------:|------:|:--------------------------------|
| Withaferin A             | CID265237             | 0.425 | 0.025 | Tubulin antagonist              |
|                          |                       | 0.65  | 0.012 | Apoptosis agonist               |
|                          |                       | 0.55  | 0.018 | Cytostatic                      |
|                          |                       | 0.28  | 0.035 | Microtubule formation inhibitor |
|                          |                       | 0.72  | 0.005 | Anticarcinogenic                |
| Withanolide A            | CID11294368           | 0.422 | 0.026 | Tubulin antagonist              |
|                          |                       | 0.58  | 0.015 | Apoptosis agonist               |
|                          |                       | 0.51  | 0.022 | Cytostatic                      |
|                          |                       | 0.68  | 0.008 | Anticarcinogenic                |
| alpha-Glycyrrhizin       | CID158471             | 0.397 | 0.045 | Apoptosis agonist               |
|                          |                       | 0.45  | 0.032 | Anticarcinogenic                |
|                          |                       | 0.35  | 0.042 | Cytostatic                      |
| 27-O-acetyl-Withaferin A | CID57328756           | 0.39  | 0.032 | Tubulin antagonist              |
|                          |                       | 0.61  | 0.018 | Apoptosis agonist               |
|                          |                       | 0.49  | 0.028 | Cytostatic                      |
|                          |                       | 0.24  | 0.045 | Microtubule formation inhibitor |
| Celastrol                | CID122724             | 0.78  | 0.004 | Apoptosis agonist               |
|                          |                       | 0.71  | 0.008 | Cytostatic                      |
|                          |                       | 0.82  | 0.002 | Anticarcinogenic                |
|                          |                       | 0.307 | 0.054 | Tubulin antagonist              |
| Tingenone                | CID101520             | 0.68  | 0.012 | Apoptosis agonist               |
|                          |                       | 0.59  | 0.019 | Cytostatic                      |
|                          |                       | 0.287 | 0.062 | Tubulin antagonist              |
| Pristimerin              | CID159516             | 0.74  | 0.006 | Apoptosis agonist               |
|                          |                       | 0.67  | 0.011 | Cytostatic                      |
|                          |                       | 0.267 | 0.071 | Tubulin antagonist              |
