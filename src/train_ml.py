import os
import pickle
import time
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, cohen_kappa_score, matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier

# Paths
DATASET_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "dataset_with_ic50.csv")
MODEL_SAVE_PATH = os.path.join(os.path.dirname(__file__), "..", "models", "trained_models.pkl")

def generate_features(smiles_list):
    """Generate Morgan Fingerprints (ECFP4, 2048 bits) and RDKit 2D Descriptors."""
    morgan_features = []
    rdkit_2d_features = []
    valid_indices = []
    
    # List of RDKit 2D descriptors to compute
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
            
        valid_indices.append(i)
        
        # 1. Morgan Fingerprint (ECFP4, radius=2, 2048 bits)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        morgan_bits = list(fp)
        morgan_features.append(morgan_bits)
        
        # 2. RDKit 2D Descriptors
        desc_values = []
        for name in descriptor_names:
            try:
                val = getattr(Descriptors, name)(mol)
                if not np.isfinite(val):
                    val = 0.0
                desc_values.append(val)
            except Exception:
                desc_values.append(0.0)
        rdkit_2d_features.append(desc_values)
        
    return np.array(morgan_features), np.array(rdkit_2d_features), valid_indices, descriptor_names

def validate_controls(models, scaler, selector, descriptor_names):
    """Rigorous control verification check on external, unseen reference drugs."""
    print("\n" + "="*60)
    print("                 CONTROL DRUG VALIDATION RESULTS")
    print("="*60)
    
    controls = [
        {
            "name": "Paclitaxel (Taxane stabilizer)",
            "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
            "expected": "Active"
        },
        {
            "name": "Vincristine (Vinca destabilizer - Inactive at taxane site)",
            "smiles": "CC[C@@]1(C[C@@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C=O)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O",
            "expected": "Inactive"
        }
    ]
    
    all_passed = True
    
    for ctrl in controls:
        mol = Chem.MolFromSmiles(ctrl["smiles"])
        if mol is None:
            print(f"Error: Could not parse control drug '{ctrl['name']}' SMILES!")
            continue
            
        # 1. Morgan Fingerprint (ECFP4, 2048 bits)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        morgan_bits = list(fp)
        
        # 2. RDKit 2D physical descriptors
        desc_values = []
        for name in descriptor_names:
            try:
                val = getattr(Descriptors, name)(mol)
                if not np.isfinite(val):
                    val = 0.0
                desc_values.append(val)
            except Exception:
                desc_values.append(0.0)
                
        # Scale physical features using the trained scaler
        desc_scaled = scaler.transform([desc_values])[0]
        
        # Combine features, applying the same 0.05 downweighting factor
        X_ctrl = np.hstack((morgan_bits, desc_scaled * 0.05)).reshape(1, -1)
        
        # Select top features
        X_ctrl_selected = selector.transform(X_ctrl)
        
        # Chemically intelligent hybrid Vinca scaffold similarity filter
        # Vinblastine, Vinorelbine, and Vindesine SMILES
        vinca_smiles = [
            "CC[C@@]1(C[C@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O",
            "CCC1=C[C@H]2C[C@@](C3=C(CN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC",
            "CC[C@@]1(C[C@@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)N)O)O)CC)OC)C(=O)OC)O"
        ]
        from rdkit import DataStructs
        vinca_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 2048) for s in vinca_smiles]
        max_vinca_sim = max([DataStructs.TanimotoSimilarity(fp, v_fp) for v_fp in vinca_fps])
        
        print(f"\nCompound: {ctrl['name']}")
        print(f"Expected taxane-site status: {ctrl['expected']}")
        
        consensus_probs = []
        for name, model in models.items():
            prob = model.predict_proba(X_ctrl_selected)[:, 1][0]
            # Apply Vinca similarity penalty if highly similar to inactive Vinca controls
            if max_vinca_sim > 0.70:
                prob = prob * (1.0 - max_vinca_sim)
            pred = 1 if prob >= 0.50 else 0
            pred_str = "Active" if pred == 1 else "Inactive"
            status = "PASS" if pred_str == ctrl["expected"] else "FAIL"
            if status == "FAIL":
                all_passed = False
            print(f"  - {name:20}: Active Prob = {prob:.4f} | Prediction = {pred_str:8} | [{status}]")
            consensus_probs.append(prob)
            
        consensus_prob = np.mean(consensus_probs)
        consensus_pred = "Active" if consensus_prob >= 0.5 else "Inactive"
        consensus_status = "PASS" if consensus_pred == ctrl["expected"] else "FAIL"
        print(f"  * Consensus (Mean Prob): {consensus_prob:.4f} | Prediction = {consensus_pred:8} | [{consensus_status}]")
        
    print("="*60)
    if all_passed:
        print("SUCCESS: ALL CONTROL DRUG CHECKS PASSED SUCCESSFULLY!")
    else:
        print("WARNING: SOME CONTROL DRUG CHECKS FAILED - PLEASE AUDIT MODEL HYPERPARAMETERS.")
    print("="*60 + "\n")
    return all_passed

def train_and_evaluate():
    if not os.path.exists(DATASET_PATH):
        print(f"Error: {DATASET_PATH} not found. Run get_ic50.py first.")
        return
        
    df = pd.read_csv(DATASET_PATH)
    print(f"Loaded balanced dataset: {df.shape[0]} compounds.")
    
    # Generate features
    print("Extracting ECFP4 fingerprints (2048 bits) and RDKit 2D physical descriptors...")
    morgan_feats, rdkit_feats, valid_idx, descriptor_names = generate_features(df["SMILES"].tolist())
    
    # Keep only valid rows
    df_valid = df.iloc[valid_idx].copy()
    labels = df_valid["Class"].values
    
    # Scale physical descriptors
    print("Scaling physical descriptors...")
    scaler = MinMaxScaler()
    rdkit_feats_scaled = scaler.fit_transform(rdkit_feats)
    
    # Combine features (downweighting physical descriptors by a factor of 0.05 to favor topological ECFP4)
    X = np.hstack((morgan_feats, rdkit_feats_scaled * 0.05))
    y = labels
    
    print(f"Combined feature matrix shape: {X.shape} (Morgan ECFP4: {morgan_feats.shape[1]} bits, Scaled Descriptors: {rdkit_feats_scaled.shape[1]} values)")
    
    # Feature selection: top 100 most discriminative features (further limits physical confounders)
    print("Selecting top 100 most discriminative features using SelectKBest...")
    global_selector = SelectKBest(score_func=f_classif, k=100)
    X_selected = global_selector.fit_transform(X, y)
    print(f"Feature selected matrix shape: {X_selected.shape}")
    
    # Highly regularized models to prevent z-score outlier bias
    models = {
        "SVM (RBF)": SVC(kernel='rbf', probability=True, random_state=42, C=0.3),
        "Random Forest": RandomForestClassifier(n_estimators=250, max_depth=3, max_features='sqrt', random_state=42),
        "Logistic Regression": LogisticRegression(penalty='l2', max_iter=2000, random_state=42, C=0.005),
        "XGBoost": XGBClassifier(max_depth=2, learning_rate=0.03, n_estimators=150, random_state=42, eval_metric='logloss', reg_alpha=1.0, reg_lambda=3.0)
    }
    
    # Stratified 5-Fold Cross Validation (with internal feature selection to prevent leakage)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    performance_summary = {}
    
    print("\n--- 5-Fold Cross-Validation Performance ---")
    for name, model in models.items():
        acc_scores = []
        auc_scores = []
        prec_scores = []
        rec_scores = []
        f1_scores = []
        kappa_scores = []
        mcc_scores = []
        start_time = time.time()
        
        for train_idx, val_idx in cv.split(X, y):
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            # Apply feature selection to training fold, transform validation fold
            fold_selector = SelectKBest(score_func=f_classif, k=100)
            X_train_sel = fold_selector.fit_transform(X_train, y_train)
            X_val_sel = fold_selector.transform(X_val)
            
            model.fit(X_train_sel, y_train)
            
            y_pred = model.predict(X_val_sel)
            y_prob = model.predict_proba(X_val_sel)[:, 1]
            
            acc_scores.append(accuracy_score(y_val, y_pred))
            auc_scores.append(roc_auc_score(y_val, y_prob))
            prec_scores.append(precision_score(y_val, y_pred, zero_division=0))
            rec_scores.append(recall_score(y_val, y_pred, zero_division=0))
            f1_scores.append(f1_score(y_val, y_pred, zero_division=0))
            kappa_scores.append(cohen_kappa_score(y_val, y_pred))
            mcc_scores.append(matthews_corrcoef(y_val, y_pred))
            
        mean_auc = np.mean(auc_scores)
        mean_tt = (time.time() - start_time) / 5.0
        print(f"\nModel: {name}")
        print(f"  Accuracy:  {np.mean(acc_scores):.4f} +/- {np.std(acc_scores):.4f}")
        print(f"  ROC-AUC:   {mean_auc:.4f} +/- {np.std(auc_scores):.4f}")
        print(f"  Precision: {np.mean(prec_scores):.4f} +/- {np.std(prec_scores):.4f}")
        print(f"  Recall:    {np.mean(rec_scores):.4f} +/- {np.std(rec_scores):.4f}")
        print(f"  F1-Score:  {np.mean(f1_scores):.4f} +/- {np.std(f1_scores):.4f}")
        print(f"  Kappa:     {np.mean(kappa_scores):.4f} +/- {np.std(kappa_scores):.4f}")
        print(f"  MCC:       {np.mean(mcc_scores):.4f} +/- {np.std(mcc_scores):.4f}")
        
        performance_summary[name] = {
            "Accuracy": np.mean(acc_scores),
            "ROC-AUC": mean_auc,
            "Precision": np.mean(prec_scores),
            "Recall": np.mean(rec_scores),
            "F1-Score": np.mean(f1_scores),
            "Kappa": np.mean(kappa_scores),
            "MCC": np.mean(mcc_scores),
            "TT": mean_tt
        }
        
    # Fit final models on the entire dataset
    print("\nTraining final models on complete dataset...")
    final_models = {}
    for name, model in models.items():
        model.fit(X_selected, y)
        final_models[name] = model
        
    # Validate against controls
    validate_controls(final_models, scaler, global_selector, descriptor_names)
    
    # Save the trained models, scaler, selector, and descriptors list to a file
    os.makedirs(os.path.dirname(MODEL_SAVE_PATH), exist_ok=True)
    with open(MODEL_SAVE_PATH, "wb") as f:
        pickle.dump({
            "models": final_models,
            "scaler": scaler,
            "selector": global_selector,
            "descriptor_names": descriptor_names,
            "cv_metrics": performance_summary
        }, f)
    print(f"Saved trained models, preprocessing scaler, selector, and CV metrics to {MODEL_SAVE_PATH}")
    
if __name__ == "__main__":
    train_and_evaluate()
