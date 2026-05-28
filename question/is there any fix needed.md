# Question: Is there any fix needed?

## Answer:

Yes. Based on the findings in the **Deep Audit Report**, there were several critical fixes required to align the code, data, and documentation with correct scientific standards and to make the pipeline robust against external service failures:

### 1. Corrected the Delaney ESOL Solubility Formula (Critical)
* **File:** `src/admet_prediction.py`
* **Issue:** The ESOL equation implementation had incorrect signs and coefficient magnitudes for *Rotatable Bonds (RB)* and *Aromatic Proportion (AP)*.
* **Fix:** Corrected the equation to align with the standard Delaney (2004) publication:
  $$\text{LogS} = 0.16 - 0.63 \cdot \text{LogP} - 0.0062 \cdot \text{MW} + 0.066 \cdot \text{RB} - 0.74 \cdot \text{AP}$$
  This updated all aqueous solubility estimations for the candidate hits and control drugs to their scientifically correct values.

### 2. Added Robust Fallback for ChEMBL API Downtime (Critical)
* **File:** `src/get_ic50.py`
* **Issue:** Public biology servers like ChEMBL occasionally encounter service downtime (returning `500 Internal Server Error`), which was causing the entire pipeline execution to abort due to a `ValueError` (sampling from an empty array).
* **Fix:** Added a robust fallback check. If the API queries return empty results but a valid, high-quality processed dataset (`dataset_with_ic50.csv` containing exactly 300 compounds) already exists locally, the script gracefully logs the event and reuses the existing dataset to prevent the pipeline from crashing.

### 3. Archived the Legacy hit-prediction Script (Major)
* **Action:** Renamed `src/predict_hits.py` to `src/predict_hits_LEGACY.py`.
* **Issue:** The legacy script was using obsolete MACCS keys and lacked the SelectKBest feature selection, physical descriptor downweighting, and Tanimoto scaffold penalties present in the current models. Running it would crash or produce corrupt predictions.
* **Fix:** Safely renamed it to avoid confusion and ensure users run the correct compound resolver, `resolve_and_predict.py`.

### 4. Synchronized Documentation & Reports (Critical)
* **Files:** `README.md` and `results/drug_discovery_report.md`
* **Issue:** Both files described a stale, pre-overhaul pipeline (e.g., a 60-compound dataset instead of 300, MACCS keys instead of 2048-bit Morgan Fingerprints, `StandardScaler` instead of `MinMaxScaler`, and claims that all 7 hits were predicted inactive).
* **Fix:** Rewrote these sections to accurately reflect the balanced 300-compound dataset size, Morgan ECFP4 fingerprints, downweighted scaling methodology, ANOVA SelectKBest, and the new prediction scores (where **5 of 7 hits are now predicted active stabilizers** by the robust ML models).
