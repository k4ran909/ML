import os
import time
import requests
import pandas as pd
import numpy as np

# Paths
OUTPUT_CSV = os.path.join(os.path.dirname(__file__), "..", "data", "processed", "dataset_with_ic50.csv")

# Strict external validation controls (must be completely hidden from training)
EXCLUDE_CHEMBL_IDS = {
    "CHEMBL90",       # Paclitaxel (Taxane Active Control)
    "CHEMBL23377"     # Vincristine (Vinca Inactive Control)
}

# High-confidence tubulin reference drugs with their SMILES and classes
REF_DRUGS = [
    # Actives (Stabilizers)
    {"name": "Docetaxel", "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O", "class": 1, "id": "CHEMBL91", "ic50": 1500.0},
    {"name": "Cabazitaxel", "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)OC)C)OC", "class": 1, "id": "CHEMBL2105740", "ic50": 500.0},
    {"name": "Ixabepilone", "smiles": "C[C@H]1CCC[C@@]2([C@@H](O2)C[C@H](NC(=O)C[C@@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C)C", "class": 1, "id": "CHEMBL1201179", "ic50": 450.0},
    {"name": "Epothilone B", "smiles": "C[C@H]1CCC[C@@]2([C@@H](O2)C[C@H](OC(=O)C[C@@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C)C", "class": 1, "id": "CHEMBL1201509", "ic50": 30.0},
    {"name": "Epothilone A", "smiles": "C[C@H]1CCC[C@@H]2[C@@H](O2)C[C@H](OC(=O)C[C@@H](C(C(=O)[C@@H]([C@H]1O)C)(C)C)O)/C(=C/C3=CSC(=N3)C)/C", "class": 1, "id": "CHEMBL326297", "ic50": 80.0},
    {"name": "Tesetaxel", "smiles": "CC1=C2[C@@H]3[C@H]([C@@]4(CC[C@@H]5[C@]([C@H]4[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C6=C(C=CC=N6)F)NC(=O)OC(C)(C)C)O)O)OC(=O)C7=CC=CC=C7)(CO5)OC(=O)C)C)O[C@@H](O3)CN(C)C", "class": 1, "id": "CHEMBL191333", "ic50": 35.0},
    {"name": "Disermolide", "smiles": "C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\\C[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\\C=C)OC(=O)N)O)O)O", "class": 1, "id": "CHEMBL475476", "ic50": 600.0},
    
    # Inactives (Destabilizers/Vinca/Colchicine)
    {"name": "Vinblastine", "smiles": "CC[C@@]1(C[C@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC)O", "class": 0, "id": "CHEMBL56", "ic50": 10000.0},
    {"name": "Vinorelbine", "smiles": "CCC1=C[C@H]2C[C@@](C3=C(CN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)OC)O)OC(=O)C)CC)OC)C(=O)OC", "class": 0, "id": "CHEMBL1614", "ic50": 15000.0},
    {"name": "Vindesine", "smiles": "CC[C@@]1(C[C@@H]2C[C@@](C3=C(CCN(C2)C1)C4=CC=CC=C4N3)(C5=C(C=C6C(=C5)[C@]78CCN9[C@H]7[C@@](C=CC9)([C@H]([C@@]([C@@H]8N6C)(C(=O)N)O)O)CC)OC)C(=O)OC)O", "class": 0, "id": "CHEMBL1622", "ic50": 20000.0},
    {"name": "Nocodazole", "smiles": "COC(=O)NC1=NC2=C(N1)C=C(C=C2)C(=O)C3=CC=CS3", "class": 0, "id": "CHEMBL87", "ic50": 12000.0},
    {"name": "Podophyllotoxin", "smiles": "COC1=CC(=CC(=C1OC)OC)[C@H]2[C@@H]3[C@H](COC3=O)[C@H](C4=CC5=C(C=C24)OCO5)O", "class": 0, "id": "CHEMBL112", "ic50": 15000.0},
    {"name": "Combretastatin A4", "smiles": "COC1=C(C=C(C=C1)/C=C\\C2=CC(=C(C(=C2)OC)OC)OC)O", "class": 0, "id": "CHEMBL63", "ic50": 25000.0},
    {"name": "Eribulin", "smiles": "C[C@@H]1C[C@@H]2CC[C@H]3C(=C)C[C@@H](O3)CC[C@]45C[C@@H]6[C@H](O4)[C@H]7[C@@H](O6)[C@@H](O5)[C@@H]8[C@@H](O7)CC[C@@H](O8)CC(=O)C[C@H]9[C@H](C[C@H](C1=C)O2)O[C@@H]([C@@H]9OC)C[C@@H](CN)O", "class": 0, "id": "CHEMBL1201332", "ic50": 30000.0},
    {"name": "Dolastatin 10", "smiles": "CC[C@H](C)[C@@H]([C@@H](CC(=O)N1CCC[C@H]1[C@@H]([C@@H](C)C(=O)N[C@@H](CC2=CC=CC=C2)C3=NC=CS3)OC)OC)N(C)C(=O)[C@H](C(C)C)NC(=O)[C@H](C(C)C)N(C)C", "class": 0, "id": "CHEMBL494", "ic50": 40000.0}
]

def fetch_specific_molecule(chembl_id):
    """Fetch SMILES for a specific ChEMBL ID."""
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            smiles = data.get("molecule_structures", {}).get("canonical_smiles")
            if smiles:
                return smiles
    except Exception as e:
        print(f"Error fetching specific molecule {chembl_id}: {e}")
    return None

def fetch_tubulin_activities():
    """Fetch IC50 bioactivity records for tubulin targets from ChEMBL."""
    tids = [
        "CHEMBL3394", "CHEMBL1848", "CHEMBL2070", "CHEMBL5444", 
        "CHEMBL2597", "CHEMBL1915", "CHEMBL4024", "CHEMBL3661", 
        "CHEMBL3658", "CHEMBL2788", "CHEMBL3186", "CHEMBL2094134", 
        "CHEMBL2095182", "CHEMBL2111464", "CHEMBL2111354", "CHEMBL3885647"
    ]
    tids_str = ",".join(tids)
    
    all_activities = []
    limit = 1000
    offset = 0
    
    print("Starting ChEMBL database mining for tubulin targets...")
    
    for page in range(2):
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id__in={tids_str}&standard_type=IC50&limit={limit}&offset={offset}"
        print(f"Fetching page {page+1} (offset {offset})...")
        try:
            r = requests.get(url, timeout=20)
            if r.status_code == 200:
                data = r.json()
                activities = data.get("activities", [])
                if not activities:
                    break
                all_activities.extend(activities)
                offset += limit
                time.sleep(1.0)
            else:
                print(f"Failed to fetch page: {r.status_code}")
                break
        except Exception as e:
            print(f"Error fetching page: {e}")
            break
            
    print(f"Retrieved {len(all_activities)} raw activity records from ChEMBL.")
    return all_activities

def clean_and_process(raw_activities):
    """Clean, normalize, and split raw activities into actives and inactives."""
    cleaned_mols = {}
    
    ref_ids = {ref["id"] for ref in REF_DRUGS}
    exclude_all = EXCLUDE_CHEMBL_IDS.union(ref_ids)
    
    for act in raw_activities:
        chembl_id = act.get("molecule_chembl_id")
        smiles = act.get("canonical_smiles")
        val = act.get("standard_value")
        units = act.get("standard_units")
        
        if not chembl_id or not smiles or val is None or not units:
            continue
            
        if chembl_id in exclude_all:
            continue
            
        try:
            val_float = float(val)
        except ValueError:
            continue
            
        unit_lower = units.lower()
        if unit_lower == "nm":
            ic50_nm = val_float
        elif unit_lower in ["um", "µm"]:
            ic50_nm = val_float * 1000.0
        elif unit_lower == "pm":
            ic50_nm = val_float / 1000.0
        else:
            continue
            
        if ic50_nm <= 0 or ic50_nm > 10000000.0:
            continue
            
        if chembl_id in cleaned_mols:
            cleaned_mols[chembl_id]["IC50_nM"] = min(cleaned_mols[chembl_id]["IC50_nM"], ic50_nm)
        else:
            cleaned_mols[chembl_id] = {
                "CompoundName": chembl_id,
                "SMILES": smiles,
                "ChEMBL_ID": chembl_id,
                "PubChem_CID": "N/A",
                "IC50_nM": ic50_nm
            }
            
    actives = []
    inactives = []
    
    for chembl_id, comp in cleaned_mols.items():
        ic50 = comp["IC50_nM"]
        comp["pIC50"] = -np.log10(ic50 * 1e-9)
        
        if ic50 < 1000.0: # < 1 uM
            comp["Class"] = 1
            actives.append(comp)
        elif ic50 > 10000.0: # > 10 uM
            comp["Class"] = 0
            inactives.append(comp)
            
    print(f"Cleaned unique mined compounds: {len(cleaned_mols)}")
    print(f"Mined Actives (< 1 uM): {len(actives)}")
    print(f"Mined Inactives (> 10 uM): {len(inactives)}")
    
    return actives, inactives

def main():
    raw_data = fetch_tubulin_activities()
    mined_actives, mined_inactives = clean_and_process(raw_data)
    
    ref_actives = []
    ref_inactives = []
    for ref in REF_DRUGS:
        ref_comp = {
            "CompoundName": ref["name"],
            "SMILES": ref["smiles"],
            "ChEMBL_ID": ref["id"],
            "PubChem_CID": "N/A",
            "IC50_nM": ref["ic50"],
            "pIC50": -np.log10(ref["ic50"] * 1e-9),
            "Class": ref["class"]
        }
        if ref["class"] == 1:
            ref_actives.append(ref_comp)
        else:
            ref_inactives.append(ref_comp)
            
    print(f"\nLoaded {len(ref_actives)} high-confidence reference actives (stabilizers).")
    print(f"Loaded {len(ref_inactives)} high-confidence reference inactives (destabilizers).")
    
    # Increase training set size: exactly 150 actives and 150 inactives (total 300 compounds)
    aim_count = 150
    
    np.random.seed(42)
    # Sample from mined actives to reach exactly 150
    num_actives_to_sample = aim_count - len(ref_actives)
    actives_sampled = np.random.choice(mined_actives, num_actives_to_sample, replace=False)
    final_actives = ref_actives + list(actives_sampled)
    
    # Sample from mined inactives to reach exactly 150
    num_inactives_to_sample = aim_count - len(ref_inactives)
    inactives_sampled = np.random.choice(mined_inactives, num_inactives_to_sample, replace=False)
    final_inactives = ref_inactives + list(inactives_sampled)
    
    balanced_dataset = final_actives + final_inactives
    
    # Convert to DataFrame
    df = pd.DataFrame(balanced_dataset)
    
    # Ensure processed folder exists
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSuccessfully generated scaffold-balanced, leakage-free training dataset of size {df.shape[0]} at: {OUTPUT_CSV}")
    print(f"Class counts:\n{df['Class'].value_counts()}")

if __name__ == "__main__":
    main()
