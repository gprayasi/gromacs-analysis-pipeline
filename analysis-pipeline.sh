#!/bin/bash

################################################################################
#                                                                              #
#             UNIFIED GROMACS MOLECULAR DYNAMICS ANALYSIS PIPELINE             #
#                                                                              #
################################################################################
#
# DESCRIPTION:
# This script automates the post-simulation analysis of GROMACS trajectories.
# It performs a series of standard and advanced analyses, including RMSD, RMSF,
# radius of gyration, secondary structure, principal component analysis (PCA),
# and free energy landscape (FEL) generation.
#
# It is designed to be reusable for different simulation systems by changing
# the variables in the 'USER CONFIGURATION' section.
#
#
# USAGE:
# 1. Make the script executable:
#    chmod +x analysis_pipeline.sh
#
# 2. Run the script for your system:
#    ./analysis_pipeline.sh
#
# 3. Ensure all required input files and external scripts are in the
#    working directory or correctly pathed.
#
#
# REQUIREMENTS:
# [cite_start]- GROMACS (tested with version 2022) [cite: 8, 18, 86]
# - DSSP (and the 'DSSP' environment variable set)
# - Perl (for the 'sham.pl' script)
# - Python 2.7 (for the 'xpm2txt.py' script)
# - Python 3 (for the 'fel3d.py' script, if used)
#
#
# REQUIRED FILES IN WORKING DIRECTORY:
# - A GROMACS TPR file (e.g., 'system.tpr')
# - A GROMACS XTC trajectory file (e.g., 'system.xtc')
# - 'sham.pl': A perl script to prepare data for gmx sham.
# - 'xpm2txt.py': A python script to convert XPM matrix to a text file.
# - 'parameters.m2p': A parameter file for gmx xpm2ps for FEL plotting (optional).
# - 'dssp.m2p': A parameter file for gmx xpm2ps for DSSP plotting (optional).
#
################################################################################

set -e # Exit immediately if a command exits with a non-zero status.

# --- USER CONFIGURATION ---

# Prefix for all output files. This will also be used to find the input files.
# Example: For inputs 'Lac-neutral.tpr' and 'Lac-neutral.xtc', set PREFIX="Lac-neutral"
PREFIX="Lac-neutral-50ns"

# Input files. The script will look for ${PREFIX}.tpr and ${PREFIX}.xtc
TPR_FILE="${PREFIX}.tpr"
XTC_FILE="${PREFIX}.xtc"

# Residue ranges for creating custom index groups.
# These are based on the examples in the command history. Modify as needed.
[cite_start]FLUCTUATION_ZONE_RESIDUES="r 360-420" # [cite: 75, 403, 620, 778, 1002]
[cite_start]DOMAIN1_RESIDUES="r 46-164" # [cite: 234]
[cite_start]DOMAIN2_RESIDUES="r 179-321" # [cite: 235]
[cite_start]DOMAIN3_RESIDUES="r 413-519" # [cite: 236]
BINDING_SITE_RESIDUES="r 101 103 141 143 446 449 451 503 504 502 508"

# Principal components to use for Free Energy Landscape.
# Inspect the cosine content from Step 5 to choose the best PCs.
PC_FOR_FEL_1=2
PC_FOR_FEL_2=3

# --- END OF USER CONFIGURATION ---


# --- CHECK FOR INPUT FILES ---
if ! [ -f "$TPR_FILE" ] || ! [ -f "$XTC_FILE" ]; then
    echo "ERROR: Input files not found!"
    echo "Please ensure '${TPR_FILE}' and '${XTC_FILE}' are in the current directory."
    exit 1
fi

echo "Pipeline started for system: $PREFIX"


################################################################################
# SECTION 1: TRAJECTORY PRE-PROCESSING
################################################################################
echo "### Section 1: Trajectory Pre-processing ###"

NO_PBC_XTC="${PREFIX}_noPBC.xtc"
if [ ! -f "$NO_PBC_XTC" ]; then
    echo "Step 1.1: Removing periodic boundary conditions and centering..."
    # Select group 1 ('Protein') for centering and group 0 ('System') for output.
    [cite_start]echo -e "1\n0" | gmx trjconv -s "$TPR_FILE" -f "$XTC_FILE" -o "$NO_PBC_XTC" -pbc mol -center -tu ns [cite: 1]
    echo "Centering complete."
else
    echo "Step 1.1: Found existing no-PBC trajectory file. Skipping."
fi


################################################################################
# SECTION 2: STANDARD ANALYSIS
################################################################################
echo -e "\n### Section 2: Standard Analysis ###"

# --- RMSD ---
echo "Step 2.1: Calculating RMSD of protein backbone..."
# Select group 4 ('Backbone') for both fitting and calculation.
[cite_start]echo -e "4\n4" | gmx rms -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_rmsd_backbone.xvg" -tu ns [cite: 8]

# --- RMSF ---
echo "Step 2.2: Calculating RMSF of C-alpha atoms..."
# Select group 3 ('C-alpha') for calculation.
[cite_start]echo -e "3" | gmx rmsf -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_rmsf_calpha.xvg" -res [cite: 15]

# --- Radius of Gyration ---
echo "Step 2.3: Calculating Radius of Gyration..."
# Select group 1 ('Protein').
[cite_start]echo -e "1" | gmx gyrate -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_gyrate.xvg" [cite: 77, 261, 393]

# --- Hydrogen Bonds ---
echo "Step 2.4: Calculating intra-protein hydrogen bonds..."
# Select group 1 ('Protein') for both groups.
[cite_start]echo -e "1\n1" | gmx hbond -s "$TPR_FILE" -f "$NO_PBC_XTC" -num "${PREFIX}_hbond_intra_protein.xvg" [cite: 184, 609]


################################################################################
# SECTION 3: SECONDARY STRUCTURE (DSSP)
################################################################################
echo -e "\n### Section 3: Secondary Structure Analysis ###"

if [ -z "$DSSP" ]; then
    export DSSP=$(which dssp)
fi

if [ -z "$DSSP" ] || ! [ -x "$DSSP" ]; then
    echo "WARNING: DSSP not found or not executable. Skipping secondary structure analysis."
else
    echo "Step 3.1: Calculating secondary structure with DSSP..."
    # Select group 1 ('Protein').
    [cite_start]echo -e "1" | gmx do_dssp -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_dssp.xpm" -tu ns -dt 1 [cite: 264, 461, 638]

    echo "Step 3.2: Converting DSSP output to EPS..."
    [cite_start]gmx xpm2ps -f "${PREFIX}_dssp.xpm" -o "${PREFIX}_dssp.eps" -rainbow red [cite: 85, 269, 467]
fi


################################################################################
# SECTION 4: PRINCIPAL COMPONENT ANALYSIS (PCA)
################################################################################
echo -e "\n### Section 4: Principal Component Analysis (PCA) ###"

# --- Covariance Analysis ---
echo "Step 4.1: Calculating covariance matrix for protein backbone..."
# Select group 4 ('Backbone') for both fitting and analysis.
[cite_start]echo -e "4\n4" | gmx covar -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_eigenval.xvg" -av "${PREFIX}_average.pdb" -xpm "${PREFIX}_covar.xpm" -v "${PREFIX}_eigenvec.trr" [cite: 24, 86, 324]

# --- Eigenvector Projections ---
echo "Step 4.2: Projecting trajectory onto first 5 eigenvectors..."
for i in {1..5}; do
    echo "Projecting onto PC${i}..."
    # Select group 4 ('Backbone') for both.
    echo -e "4\n4" | gmx anaeig -s "$TPR_FILE" -f "$NO_PBC_XTC" -v "${PREFIX}_eigenvec.trr" -eig "${PREFIX}_eigenval.xvg" \
    [cite_start]-proj "${PREFIX}_proj_pc${i}.xvg" -extr "${PREFIX}_extreme_pc${i}.pdb" -first "$i" -last "$i" -nframes 100 [cite: 34, 45, 58]
done

# --- Cosine Content Analysis ---
echo "Step 4.3: Analyzing cosine content of projections (to check for convergence)..."
for i in {1..5}; do
    [cite_start]gmx analyze -f "${PREFIX}_proj_pc${i}.xvg" -cc "${PREFIX}_proj_pc${i}_cc.xvg" > "${PREFIX}_proj_pc${i}_cosine.txt" 2>&1 [cite: 43, 97, 109]
done
echo "Cosine content analysis complete. Check *_cosine.txt files."
echo "Low cosine values (<0.2) suggest the projection has converged."


################################################################################
# SECTION 5: FREE ENERGY LANDSCAPE (FEL)
################################################################################
echo -e "\n### Section 5: Free Energy Landscape (FEL) ###"

# --- 2D Projection ---
echo "Step 5.1: Generating 2D projection for PC${PC_FOR_FEL_1} and PC${PC_FOR_FEL_2}..."
# Select group 4 ('Backbone') for both.
echo -e "4\n4" | gmx anaeig -s "$TPR_FILE" -f "$NO_PBC_XTC" -v "${PREFIX}_eigenvec.trr" -eig "${PREFIX}_eigenval.xvg" \
[cite_start]-2d "${PREFIX}_proj_scatter_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.xvg" -first "${PC_FOR_FEL_1}" -last "${PC_FOR_FEL_2}" [cite: 166, 372]

# --- Gibbs Free Energy Landscape ---
echo "Step 5.2: Calculating Gibbs free energy landscape with gmx sham..."
[cite_start]gmx sham -f "${PREFIX}_proj_scatter_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.xvg" -ls "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.xpm" [cite: 175]

# --- FEL Post-processing ---
echo "Step 5.3: Post-processing FEL data..."
if ! [ -x "$(command -v python2.7)" ] || ! [ -f "xpm2txt.py" ]; then
    echo "WARNING: python2.7 or xpm2txt.py not found. Skipping FEL text conversion."
else
    python2.7 xpm2txt.py -f "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.xpm" -o "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.txt"
    sort -k3,3n "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.txt" > "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}_sorted.txt"
fi
[cite_start]gmx xpm2ps -f "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.xpm" -o "${PREFIX}_fel_pc${PC_FOR_FEL_1}-${PC_FOR_FEL_2}.eps" -rainbow red [cite: 176, 381]


################################################################################
# SECTION 6: OPTIONAL / SPECIALIZED ANALYSIS
# (Uncomment the blocks you wish to run)
################################################################################
echo -e "\n### Section 6: Optional Analyses (currently commented out) ###"

# --- Analysis of Specific Protein Domains ---
# Creates an index file with three domains and calculates RMSD for each.
# echo "Step 6.1: Analyzing protein domain RMSD..."
# [cite_start]echo -e "${DOMAIN1_RESIDUES}\nname 16 dom1\n${DOMAIN2_RESIDUES}\nname 17 dom2\n${DOMAIN3_RESIDUES}\nname 18 dom3\nq" | gmx make_ndx -f "$TPR_FILE" -o "${PREFIX}_domains.ndx" [cite: 229]
#
# [cite_start]echo "Backbone dom1" | gmx rms -s "$TPR_FILE" -f "$NO_PBC_XTC" -n "${PREFIX}_domains.ndx" -o "${PREFIX}_rmsd_dom1.xvg" -tu ns [cite: 237]
# [cite_start]echo "Backbone dom2" | gmx rms -s "$TPR_FILE" -f "$NO_PBC_XTC" -n "${PREFIX}_domains.ndx" -o "${PREFIX}_rmsd_dom2.xvg" -tu ns [cite: 245]
# [cite_start]echo "Backbone dom3" | gmx rms -s "$TPR_FILE" -f "$NO_PBC_XTC" -n "${PREFIX}_domains.ndx" -o "${PREFIX}_rmsd_dom3.xvg" -tu ns [cite: 253]
# echo "Domain analysis complete."


# --- Analysis of a High-Fluctuation Region ---
# Creates an index for a specific region and calculates its RMSD.
# echo "Step 6.2: Analyzing a high-fluctuation zone..."
# [cite_start]echo -e "${FLUCTUATION_ZONE_RESIDUES}\nname 16 fluctuation_zone\nq" | gmx make_ndx -f "$TPR_FILE" -o "${PREFIX}_fluctuation_zone.ndx" [cite: 70, 397]
#
# Fit to backbone, calculate RMSD for the fluctuation zone.
# [cite_start]echo "Backbone fluctuation_zone" | gmx rms -s "$TPR_FILE" -f "$NO_PBC_XTC" -n "${PREFIX}_fluctuation_zone.ndx" -o "${PREFIX}_rmsd_fluctuation_zone.xvg" -tu ns [cite: 222, 405]
# echo "Fluctuation zone analysis complete."


# --- Analysis of Water Molecules near a Binding Site ---
# A multi-step analysis to find and characterize water molecules in a binding site.
# echo "Step 6.3: Analyzing water molecules in the binding site..."
# # 1. Create an index file for the binding site residues.
# [cite_start]echo -e "${BINDING_SITE_RESIDUES}\nname 16 binding_resi\nq" | gmx make_ndx -f "$TPR_FILE" -o "${PREFIX}_binding_resi.ndx" [cite: 189]
#
# # 2. Order the trajectory to keep water molecules whole.
# [cite_start]echo "binding_resi TIP3" | gmx trjorder -s "$TPR_FILE" -f "$NO_PBC_XTC" -n "${PREFIX}_binding_resi.ndx" -o "${PREFIX}_ordered.xtc" -da 3.0 -r 5.0 [cite: 195]
#
# # 3. Select water molecules that are within 0.5 nm of the binding site throughout the trajectory.
# [cite_start]echo 'group "TIP3" and within 0.5 of group "binding_resi"' | gmx select -s "$TPR_FILE" -f "${PREFIX}_ordered.xtc" -n "${PREFIX}_binding_resi.ndx" -on "${PREFIX}_first_shell.ndx" [cite: 210]
#
# # 4. Calculate hydrogen bonds between the binding site and surrounding water.
# [cite_start]echo "binding_resi TIP3" | gmx hbond -s "$TPR_FILE" -f "$NO_PBC_XTC" -n "${PREFIX}_binding_resi.ndx" -num "${PREFIX}_hbond_bs-water.xvg" [cite: 216]
# echo "Binding site water analysis complete."


################################################################################
# SECTION 7: STRUCTURE & TRAJECTORY EXTRACTION
################################################################################
echo -e "\n### Section 7: Structure Extraction ###"

echo "Step 7.1: Extracting PDB structures at various time points..."
# Select group 1 ('Protein') for output.
[cite_start]echo "1" | gmx trjconv -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_0ns.pdb" -dump 0 [cite: 157, 565]
[cite_start]echo "1" | gmx trjconv -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_10ns.pdb" -dump 10000 [cite: 152, 548]
[cite_start]echo "1" | gmx trjconv -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_20ns.pdb" -dump 20000 [cite: 553, 800]
[cite_start]echo "1" | gmx trjconv -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_30ns.pdb" -dump 30000 [cite: 421, 557]
[cite_start]echo "1" | gmx trjconv -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_40ns.pdb" -dump 40000 [cite: 412, 561]
[cite_start]echo "1" | gmx trjconv -s "$TPR_FILE" -f "$NO_PBC_XTC" -o "${PREFIX}_50ns.pdb" -dump 50000 [cite: 147, 417]

echo -e "\n################################################"
echo "###      GROMACS Analysis Pipeline Finished      ###"
echo "################################################"
echo "All output files are prefixed with: $PREFIX"
