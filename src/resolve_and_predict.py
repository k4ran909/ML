import os
import urllib.parse
import requests
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

MODEL_SAVE_PATH = os.path.join(os.path.dirname(__file__), "..", "models", "trained_models.pkl")
OUTPUT_PREDICTIONS_CSV = os.path.join(os.path.dirname(__file__), "..", "results", "hit_predictions_results.csv")

# The hits provided by the user
HIT_NAMES = [
    "Withanolide A",
    "27-O-acetyl-Withaferin A",
    "Pristimerin",
    "Celastrol",
    "Withaferin A",
    "alpha-Glycyrrhizin",
    "Tingenone"
]

# Alternate search names for PubChem lookup in case the primary name fails
ALT_NAMES = {
    "27-O-acetyl-Withaferin A": ["27-O-acetylwithaferin A", "Acetylwithaferin A", "27-acetylwithaferin A"],
    "alpha-Glycyrrhizin": ["Glycyrrhizin", "Glycyrrhizic acid", "alpha-Glycyrrhizic acid", "Glycyron"],
}

def resolve_name_to_smiles(name):
    """Fetch Canonical, Isomeric, or Connectivity SMILES from PubChem for a compound name."""
    names_to_try = [name]
    if name in ALT_NAMES:
        names_to_try.extend(ALT_NAMES[name])
        
    for query_name in names_to_try:
        print(f"Trying to resolve '{query_name}' on PubChem...")
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(query_name)}/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                properties = data.get("PropertyTable", {}).get("Properties", [])
                if properties:
                    props = properties[0]
                    smiles = props.get("CanonicalSMILES") or props.get("IsomericSMILES") or props.get("ConnectivitySMILES")
                    cid = props.get("CID")
                    if smiles:
                        print(f"Success! Resolved '{query_name}' to CID: {cid}")
                        return smiles, cid
        except Exception as e:
            print(f"Network error resolving '{query_name}': {e}")
            
    return None, None

def generate_features_for_hits(smiles_list, descriptor_names):
    morgan_features = []
    rdkit_2d_features = []
    valid_indices = []
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        valid_indices.append(i)
        
        # 1. Morgan Fingerprint (ECFP4, 2048 bits)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        morgan_bits = list(fp)
        morgan_features.append(morgan_bits)
        
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
        rdkit_2d_features.append(desc_values)
        
    return np.array(morgan_features), np.array(rdkit_2d_features), valid_indices

def main():
    if not os.path.exists(MODEL_SAVE_PATH):
        print(f"Error: {MODEL_SAVE_PATH} not found. Run train_ml.py first.")
        return
        
    # Load model package
    with open(MODEL_SAVE_PATH, "rb") as f:
        saved_data = pickle.load(f)
        
    models = saved_data["models"]
    scaler = saved_data["scaler"]
    selector = saved_data["selector"]
    descriptor_names = saved_data["descriptor_names"]
    
    resolved_hits = []
    
    # Resolve all hits to SMILES
    for name in HIT_NAMES:
        smiles, cid = resolve_name_to_smiles(name)
        if smiles:
            resolved_hits.append({
                "name": name,
                "smiles": smiles,
                "cid": cid
            })
        else:
            print(f"Error: Could not resolve '{name}' to a chemical structure on PubChem.")
            
    if not resolved_hits:
        print("Error: No compounds could be resolved. Exiting.")
        return
        
    smiles_list = [h["smiles"] for h in resolved_hits]
    names_list = [h["name"] for h in resolved_hits]
    cids_list = [h["cid"] for h in resolved_hits]
    
    # Extract features
    morgan_feats, rdkit_feats, valid_idx = generate_features_for_hits(smiles_list, descriptor_names)
    
    # Scale physical features using the saved scaler
    rdkit_feats_scaled = scaler.transform(rdkit_feats)
    
    # Combine features applying same factor (0.05)
    X = np.hstack((morgan_feats, rdkit_feats_scaled * 0.05))
    
    # Select top features using the saved selector
    X_selected = selector.transform(X)
    
    # Vinca inactives for scaffold similarity check
    vinca_smiles = [
        "CC[C@@]1(C[C@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O",
        "CCC1=C[C@H]2C[C@@](C3=C(CN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC",
        "CC[C@@]1(C[C@@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)N)O)O)CC)OC)C(=O)OC)O"
    ]
    from rdkit import DataStructs
    vinca_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, 2048) for s in vinca_smiles]
    
    # Predict
    predictions = []
    for idx, orig_idx in enumerate(valid_idx):
        comp_name = names_list[orig_idx]
        smiles = smiles_list[orig_idx]
        cid = cids_list[orig_idx]
        
        # Calculate maximum similarity to Vinca inactives
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        max_vinca_sim = max([DataStructs.TanimotoSimilarity(fp, v_fp) for v_fp in vinca_fps])
        
        comp_results = {
            "CompoundName": comp_name,
            "PubChem_CID": cid,
            "SMILES": smiles
        }
        
        probs = []
        for model_name, model in models.items():
            prob = model.predict_proba(X_selected[idx:idx+1])[:, 1][0]
            # Apply Vinca similarity penalty if highly similar to Vinca inactives
            if max_vinca_sim > 0.70:
                prob = prob * (1.0 - max_vinca_sim)
            pred = 1 if prob >= 0.50 else 0
            comp_results[f"{model_name}_Prob"] = prob
            comp_results[f"{model_name}_Prediction"] = "Active" if pred == 1 else "Inactive"
            probs.append(prob)
            
        comp_results["Consensus_Active_Prob"] = np.mean(probs)
        # Final classification based on consensus (threshold 0.5)
        comp_results["Consensus_Class"] = "Active" if comp_results["Consensus_Active_Prob"] >= 0.5 else "Inactive"
        predictions.append(comp_results)
        
    df = pd.DataFrame(predictions)
    # Sort by consensus active probability descending
    df = df.sort_values(by="Consensus_Active_Prob", ascending=False)
    
    # Save predictions to CSV
    os.makedirs(os.path.dirname(OUTPUT_PREDICTIONS_CSV), exist_ok=True)
    df.to_csv(OUTPUT_PREDICTIONS_CSV, index=False)
    print(f"\nFinal prediction results saved to: {OUTPUT_PREDICTIONS_CSV}")
    
    # Print formatted output
    print("\n" + "="*90)
    print("                      ML HIT SCREENING REPORT ( Morgan ECFP4 + RDKit 2D )")
    print("="*90)
    
    cols = ["CompoundName", "PubChem_CID", "Consensus_Active_Prob", "Consensus_Class"]
    print(df[cols].to_string(index=False))
    print("="*90)

if __name__ == "__main__":
    main()
