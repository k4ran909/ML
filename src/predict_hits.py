import os
import sys
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors

# Paths
MODEL_SAVE_PATH = os.path.join(os.path.dirname(__file__), "..", "models", "trained_models.pkl")

def generate_features_for_hits(smiles_list, descriptor_names):
    """Generate MACCS keys and RDKit 2D Descriptors for screening hits."""
    maccs_features = []
    rdkit_2d_features = []
    valid_indices = []
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Invalid SMILES ignored at index {i}: {smiles}")
            continue
            
        valid_indices.append(i)
        
        # 1. MACCS Keys
        maccs = list(MACCSkeys.GenMACCSKeys(mol))[1:]
        maccs_features.append(maccs)
        
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
        
    return np.array(maccs_features), np.array(rdkit_2d_features), valid_indices

def predict(hits_data):
    """
    Predict activity of hits.
    hits_data: List of dicts with 'name' and 'smiles' keys.
    """
    if not os.path.exists(MODEL_SAVE_PATH):
        print(f"Error: {MODEL_SAVE_PATH} not found. Run train_ml.py first.")
        return
        
    with open(MODEL_SAVE_PATH, "rb") as f:
        saved_data = pickle.load(f)
        
    models = saved_data["models"]
    scaler = saved_data["scaler"]
    descriptor_names = saved_data["descriptor_names"]
    
    names = [h["name"] for h in hits_data]
    smiles_list = [h["smiles"] for h in hits_data]
    
    # Generate features
    maccs_feats, rdkit_feats, valid_idx = generate_features_for_hits(smiles_list, descriptor_names)
    
    if len(valid_idx) == 0:
        print("Error: No valid SMILES found in the hits list.")
        return
        
    # Scale physical features using the saved scaler
    rdkit_feats_scaled = scaler.transform(rdkit_feats)
    
    # Combine features
    X = np.hstack((maccs_feats, rdkit_feats_scaled))
    
    # Run prediction for each model
    results = []
    for idx, orig_idx in enumerate(valid_idx):
        comp_name = names[orig_idx]
        smiles = smiles_list[orig_idx]
        
        comp_results = {
            "CompoundName": comp_name,
            "SMILES": smiles
        }
        
        consensus_probs = []
        
        for name, model in models.items():
            prob = model.predict_proba(X[idx:idx+1])[:, 1][0]
            pred = model.predict(X[idx:idx+1])[0]
            comp_results[f"{name}_Prob"] = prob
            comp_results[f"{name}_Class"] = "Active" if pred == 1 else "Inactive"
            consensus_probs.append(prob)
            
        comp_results["Consensus_Active_Prob"] = np.mean(consensus_probs)
        results.append(comp_results)
        
    df = pd.DataFrame(results)
    # Sort by consensus active probability descending
    df = df.sort_values(by="Consensus_Active_Prob", ascending=False)
    
    # Save predictions
    output_path = os.path.join(os.path.dirname(__file__), "..", "results", "hits_ml_predictions.csv")
    df.to_csv(output_path, index=False)
    print(f"\nPredictions saved to {output_path}")
    
    # Print formatted table
    print("\n--- ML Screening Predictions Report ---")
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    print(df[["CompoundName", "Consensus_Active_Prob"] + [f"{m}_Class" for m in models.keys()]])

if __name__ == "__main__":
    # Example test usage if run directly
    test_hits = [
        {"name": "ZINC03847075", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}, # Aspirin structure for testing
        {"name": "ZINC12889138", "smiles": "CN1CCC2=CC3=C(C(=C2[C@@H]1[C@@H]4C5=C(C(=C(C=C5)OC)OC)C(=O)O4)OC)OCO3"} # Noscapine for testing
    ]
    predict(test_hits)
