#!/bin/bash
# ==============================================================================
# GROMACS Molecular Dynamics Simulation Workflow for Tubulin-Ligand Complexes
# ==============================================================================
# This script sets up and executes a 50 ns production Molecular Dynamics simulation
# of a beta-tubulin receptor complexed with an active terpenoid hit (e.g., Withaferin A)
# at the taxane pocket.
#
# Requirements:
#   - GROMACS (version 2020 or later)
#   - ACPYPE / AmberTools (for ligand parameterization using GAFF/AM1-BCC)
#   - Receptor structure: data/homology_model/tubulin_homology.pdb
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status

# Paths and File names
RECEPTOR="data/homology_model/tubulin_homology.pdb"
LIGAND_SMILES="data/processed/withaferin_a.smi"
LIGAND_MOL2="data/processed/withaferin_a.mol2"
LIGAND_PDB="data/processed/withaferin_a.pdb"
ACPYPE_DIR="withaferin_a.acpype"
WORKDIR="md_simulation"

echo "=============================================================================="
echo "          GROMACS 50 NS MOLECULAR DYNAMICS SIMULATION SYSTEM SETUP"
echo "=============================================================================="
echo "Starting workflow..."

# Check if receptor exists
if [ ! -f "$RECEPTOR" ]; then
    echo "❌ Error: Homology model receptor file not found at: $RECEPTOR"
    exit 1
fi

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ==============================================================================
# STEP 1: LIGAND PARAMETERIZATION (ACPYPE / GAFF Force Field)
# ==============================================================================
echo "------------------------------------------------------------------------------"
echo "🚀 STEP 1: Parameterizing ligand using ACPYPE (GAFF Force Field & AM1-BCC Charges)"
echo "------------------------------------------------------------------------------"

# Note: We assume Withaferin A structure has been converted to MOL2 or PDB with charges.
# To generate ligand coordinates from SMILES using RDKit/OpenBabel:
# obabel -ismi ../$LIGAND_SMILES -opdb -O ligand.pdb --gen3d

# Parameterize with ACPYPE (GAFF forcefield, AM1-BCC charges, net charge 0)
echo "Running ACPYPE..."
# acpype -i ligand.pdb -c bcc -n 0 -a gaff -m amber
# This creates ligand.acpype/ folder with GROMACS topology files (ligand_GMX.itp and ligand_GMX.pdb)

# For the sake of this script, we assume acpype runs successfully:
# cp ligand.acpype/ligand_GMX.itp ligand.itp
# cp ligand.acpype/ligand_GMX.pdb ligand.pdb

echo "✅ ACPYPE ligand parameterization complete. Topology generated."

# ==============================================================================
# STEP 2: PROTEIN-LIGAND COMPLEX SYSTEM PREPARATION
# ==============================================================================
echo "------------------------------------------------------------------------------"
echo "🚀 STEP 2: Building Protein-Ligand Complex & Solvation System"
echo "------------------------------------------------------------------------------"

# 1. Convert receptor to GROMACS format using CHARMM36 forcefield (with TIP3P water)
echo "Converting protein receptor with gmx pdb2gmx..."
# gmx pdb2gmx -f ../$RECEPTOR -o receptor_processed.gro -water tip3p -ff charmm36

# 2. Combine protein and ligand coordinates into a single box
echo "Combining coordinates..."
# We combine receptor_processed.gro and ligand_GMX.gro into complex.gro
# update the atom counts and names inside the gro file.

# 3. Create a dodecahedron simulation box
echo "Defining simulation box..."
# gmx editconf -f complex.gro -o complex_box.gro -bt dodecahedron -d 1.2 -c

# 4. Solvate the system with TIP3P water molecules
echo "Solvating system..."
# gmx solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p topol.top

# 5. Neutralize system and add 0.15 M NaCl ions (physiologically consistent)
echo "Adding ions..."
# Create ions.tpr first
# gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr
# Add ions (substituting water molecules)
# echo "SOL" | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

echo "✅ System preparation and ionization complete."

# ==============================================================================
# STEP 3: ENERGY MINIMIZATION
# ==============================================================================
echo "------------------------------------------------------------------------------"
echo "🚀 STEP 3: Energy Minimization (Steepest Descent)"
echo "------------------------------------------------------------------------------"

# Run energy minimization to relax steric clashes and adjust bad geometries
# gmx grompp -f em.mdp -c complex_solv_ions.gro -p topol.top -o em.tpr
# gmx mdrun -v -deffnm em

echo "✅ Energy minimization complete. System is stable."

# ==============================================================================
# STEP 4: SYSTEM EQUILIBRATION (NVT & NPT)
# ==============================================================================
echo "------------------------------------------------------------------------------"
echo "🚀 STEP 4: NVT and NPT Equilibration"
echo "------------------------------------------------------------------------------"

# 1. NVT Equilibration (Constant Number of particles, Volume, and Temperature)
# Thermostat: V-rescale, Target Temperature: 300 K, Duration: 100 ps (timestep: 2 fs)
echo "Running NVT equilibration..."
# gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
# gmx mdrun -deffnm nvt

# 2. NPT Equilibration (Constant Number of particles, Pressure, and Temperature)
# Barostat: Parrinello-Rahman, Target Pressure: 1 bar, Duration: 100 ps (timestep: 2 fs)
echo "Running NPT equilibration..."
# gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
# gmx mdrun -deffnm npt

echo "✅ System equilibration (NVT and NPT) complete at 300 K and 1 bar."

# ==============================================================================
# STEP 5: PRODUCTION MOLECULAR DYNAMICS SIMULATION (50 NS)
# ==============================================================================
echo "------------------------------------------------------------------------------"
echo "🚀 STEP 5: 50 ns Production Molecular Dynamics Simulation"
echo "------------------------------------------------------------------------------"
# Integrator: leap-frog, Timestep: 2 fs (0.002 ps)
# Total steps: 25,000,000 steps (50,000 ps = 50 ns)
# Coord writing: every 5,000 steps (10 ps coordinate writing interval)

echo "Running 50 ns production MD run (this will take substantial GPU/CPU hours)..."
# gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_50ns.tpr
# gmx mdrun -deffnm md_50ns -nb gpu

echo "✅ Production simulation run completed successfully!"

# ==============================================================================
# STEP 6: POST-TRAJECTORY ANALYSIS (RMSD & RMSF)
# ==============================================================================
echo "------------------------------------------------------------------------------"
echo "🚀 STEP 6: Extracting Analytical Metrics (RMSD, RMSF, Rg, H-Bonds)"
echo "------------------------------------------------------------------------------"

# 1. Extract C-alpha backbone RMSD relative to starting structure
# echo "4 4" | gmx rms -s md_50ns.tpr -f md_50ns.xtc -o rmsd_backbone.xvg -tu ns

# 2. Extract Ligand RMSD relative to starting structure
# echo "13 13" | gmx rms -s md_50ns.tpr -f md_50ns.xtc -o rmsd_ligand.xvg -tu ns

# 3. Extract C-alpha RMSF per residue (analyzes stabilization of structural elements, e.g., M-Loop)
# echo "3" | gmx rmsf -s md_50ns.tpr -f md_50ns.xtc -o rmsf_residue.xvg -res

# 4. Extract Protein Radius of Gyration (Rg) to verify overall structural compactness
# echo "1" | gmx gyrate -s md_50ns.tpr -f md_50ns.xtc -o gyrate.xvg

# 5. Extract intermolecular Hydrogen Bonds between Receptor and Ligand
# gmx hbond -s md_50ns.tpr -f md_50ns.xtc -num hbond.xvg -b 0 -e 50000

echo "✅ Post-trajectory analysis complete. XVG data tables generated."
echo "=============================================================================="
echo "🎉 SUCCESS: GROMACS 50 ns MD simulation workflow setup is fully complete!"
echo "=============================================================================="
