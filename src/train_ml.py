import os
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier

# Paths
DATASET_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "dataset_with_ic50.csv")
MODEL_SAVE_PATH = os.path.join(os.path.dirname(__file__), "..", "models", "trained_models.pkl")

def generate_features(smiles_list):
    """Generate MACCS keys and RDKit 2D Descriptors for a list of SMILES."""
    maccs_features = []
    rdkit_2d_features = []
    valid_indices = []
    
    # List of RDKit 2D descriptors to compute
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Invalid SMILES ignored at index {i}: {smiles}")
            continue
            
        valid_indices.append(i)
        
        # 1. MACCS Keys (167 bits, we ignore bit 0 as it is placeholder)
        maccs = list(MACCSkeys.GenMACCSKeys(mol))[1:]
        maccs_features.append(maccs)
        
        # 2. RDKit 2D Descriptors
        desc_values = []
        for name in descriptor_names:
            try:
                val = getattr(Descriptors, name)(mol)
                # Handle potential NaN or Inf values
                if not np.isfinite(val):
                    val = 0.0
                desc_values.append(val)
            except Exception:
                desc_values.append(0.0)
        rdkit_2d_features.append(desc_values)
        
    return np.array(maccs_features), np.array(rdkit_2d_features), valid_indices, descriptor_names

def train_and_evaluate():
    if not os.path.exists(DATASET_PATH):
        print(f"Error: {DATASET_PATH} not found. Run get_ic50.py first.")
        return
        
    df = pd.read_csv(DATASET_PATH)
    print(f"Loaded dataset: {df.shape[0]} compounds.")
    
    # Generate features
    print("Generating molecular features (MACCS Keys & RDKit 2D Descriptors)...")
    maccs_feats, rdkit_feats, valid_idx, descriptor_names = generate_features(df["SMILES"].tolist())
    
    # Keep only valid rows
    df_valid = df.iloc[valid_idx].copy()
    labels = df_valid["Class"].values
    
    # Combine features
    # Note: We scale the 2D descriptors but not the binary MACCS keys
    print("Scaling physical descriptors...")
    scaler = StandardScaler()
    rdkit_feats_scaled = scaler.fit_transform(rdkit_feats)
    
    X = np.hstack((maccs_feats, rdkit_feats_scaled))
    y = labels
    
    print(f"Final feature matrix shape: {X.shape} (MACCS: {maccs_feats.shape[1]} bits, Descriptors: {rdkit_feats_scaled.shape[1]} physical values)")
    
    # Models to train
    models = {
        "SVM (RBF)": SVC(kernel='rbf', probability=True, random_state=42),
        "Random Forest": RandomForestClassifier(n_estimators=100, max_depth=6, random_state=42),
        "Logistic Regression": LogisticRegression(penalty='l2', max_iter=1000, random_state=42),
        "XGBoost": XGBClassifier(max_depth=3, learning_rate=0.05, n_estimators=100, random_state=42, eval_metric='logloss')
    }
    
    # Stratified 5-Fold Cross Validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    performance_summary = {}
    
    print("\n--- 5-Fold Cross-Validation Performance ---")
    for name, model in models.items():
        acc_scores = []
        auc_scores = []
        prec_scores = []
        rec_scores = []
        f1_scores = []
        
        for train_idx, val_idx in cv.split(X, y):
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            model.fit(X_train, y_train)
            
            y_pred = model.predict(X_val)
            y_prob = model.predict_proba(X_val)[:, 1]
            
            acc_scores.append(accuracy_score(y_val, y_pred))
            auc_scores.append(roc_auc_score(y_val, y_prob))
            prec_scores.append(precision_score(y_val, y_pred, zero_division=0))
            rec_scores.append(recall_score(y_val, y_pred, zero_division=0))
            f1_scores.append(f1_score(y_val, y_pred, zero_division=0))
            
        print(f"\nModel: {name}")
        print(f"  Accuracy:  {np.mean(acc_scores):.4f} +/- {np.std(acc_scores):.4f}")
        print(f"  ROC-AUC:   {np.mean(auc_scores):.4f} +/- {np.std(auc_scores):.4f}")
        print(f"  Precision: {np.mean(prec_scores):.4f} +/- {np.std(prec_scores):.4f}")
        print(f"  Recall:    {np.mean(rec_scores):.4f} +/- {np.std(rec_scores):.4f}")
        print(f"  F1-Score:  {np.mean(f1_scores):.4f} +/- {np.std(f1_scores):.4f}")
        
        performance_summary[name] = {
            "Accuracy": np.mean(acc_scores),
            "ROC-AUC": np.mean(auc_scores),
            "Precision": np.mean(prec_scores),
            "Recall": np.mean(rec_scores),
            "F1-Score": np.mean(f1_scores)
        }
        
    # Fit final models on the entire dataset and save them
    print("\nTraining final models on the complete dataset...")
    final_models = {}
    for name, model in models.items():
        model.fit(X, y)
        final_models[name] = model
        
    # Save the trained models, scaler, and descriptors list to a file
    with open(MODEL_SAVE_PATH, "wb") as f:
        pickle.dump({
            "models": final_models,
            "scaler": scaler,
            "descriptor_names": descriptor_names
        }, f)
    print(f"Saved trained models and preprocessing scaler to {MODEL_SAVE_PATH}")
    
if __name__ == "__main__":
    train_and_evaluate()
