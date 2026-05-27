import os
import urllib.parse
import requests
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys
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
    maccs_features = []
    rdkit_2d_features = []
    valid_indices = []
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
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

def main():
    # Load model package
    with open(MODEL_SAVE_PATH, "rb") as f:
        saved_data = pickle.load(f)
        
    models = saved_data["models"]
    scaler = saved_data["scaler"]
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
    maccs_feats, rdkit_feats, valid_idx = generate_features_for_hits(smiles_list, descriptor_names)
    
    # Scale physical features using the saved scaler
    rdkit_feats_scaled = scaler.transform(rdkit_feats)
    
    # Combine features
    X = np.hstack((maccs_feats, rdkit_feats_scaled))
    
    # Predict
    predictions = []
    for idx, orig_idx in enumerate(valid_idx):
        comp_name = names_list[orig_idx]
        smiles = smiles_list[orig_idx]
        cid = cids_list[orig_idx]
        
        comp_results = {
            "CompoundName": comp_name,
            "PubChem_CID": cid,
            "SMILES": smiles
        }
        
        probs = []
        for model_name, model in models.items():
            prob = model.predict_proba(X[idx:idx+1])[:, 1][0]
            pred = model.predict(X[idx:idx+1])[0]
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
    df.to_csv(OUTPUT_PREDICTIONS_CSV, index=False)
    print(f"\nFinal prediction results saved to: {OUTPUT_PREDICTIONS_CSV}")
    
    # Print formatted output
    print("\n" + "="*80)
    print("                      ML HIT SCREENING REPORT")
    print("="*80)
    
    cols = ["CompoundName", "PubChem_CID", "Consensus_Active_Prob", "Consensus_Class"]
    print(df[cols].to_string(index=False))
    print("="*80)

if __name__ == "__main__":
    main()
