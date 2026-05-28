import os
import sys
import subprocess

# List of scripts to run in order of dependency
scripts = [
    ("get_ic50.py", "1. Mining experimental IC50 values from ChEMBL/PubChem"),
    ("train_ml.py", "2. Building features & training classifiers (Random Forest, SVM, LR, XGBoost)"),
    ("resolve_and_predict.py", "3. Resolving hit structures & predicting pocket binding probabilities"),
    ("admet_prediction.py", "4. Computing Lipinski properties & ESOL solubility profiling"),
    ("feature_importance.py", "5. Generating descriptor feature weights & saving importance plots"),
    ("analyze_md_trajectory.py", "6. Simulating Molecular Dynamics trajectories & computing stability RMSD/RMSF"),
    ("generate_tables.py", "7. Formatting cross-validation results & creating markdown reports"),
    ("generate_dashboard.py", "8. Building interactive HTML browser dashboard")
]

src_dir = os.path.dirname(os.path.abspath(__file__))

def main():
    print("="*75)
    print("           TUBULIN DRUG DISCOVERY PIPELINE: MASTER EXECUTION RUNNER")
    print("="*75)
    print(f"Current Interpreter: {sys.executable}")
    print(f"Working Directory: {src_dir}\n")
    
    for idx, (script, desc) in enumerate(scripts, 1):
        script_path = os.path.join(src_dir, script)
        
        if not os.path.exists(script_path):
            print(f"[ERROR] Script not found: {script_path}")
            sys.exit(1)
            
        print("="*75)
        print(f"[STEP {idx}/{len(scripts)}] Running: {desc}")
        print(f"   Command: {sys.executable} {script}")
        print("-"*75)
        
        # We execute using the exact same python interpreter (sys.executable)
        # to ensure virtual environment and package dependencies are shared.
        process = subprocess.Popen(
            [sys.executable, script],
            cwd=src_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )
        
        # Read and print the console output in real-time
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
                
        rc = process.poll()
        if rc != 0:
            print("-"*75)
            print(f"[ERROR] Step {idx} ({script}) failed with exit code {rc}.")
            print("Aborting pipeline execution to prevent downstream data corruption.")
            print("="*75)
            sys.exit(rc)
            
        print("-"*75)
        print(f"[SUCCESS] Step {idx} ({script}) completed successfully.\n")
        
    print("="*75)
    print("[SUCCESS] Entire machine learning and ADMET pipeline completed successfully!")
    print("   All results, tables, and charts are updated in the 'results/' folder.")
    print("="*75)
    print("\n" + " "*22 + "DEVELOPED & RESEARCHED BY:")
    print(" "*17 + "  ____  ___  ____ _   _    _")
    print(" "*17 + " |  _ \\|_ _/ ___| | | |  / \\")
    print(" "*17 + " | |_) || | |   | |_| | / _ \\")
    print(" "*17 + " |  _ < | | |___|  _  |/ ___ \\")
    print(" "*17 + " |_| \\_\\___\\____|_| |_/_/   \\_\\")
    print("\n" + "="*75)

if __name__ == "__main__":
    main()
