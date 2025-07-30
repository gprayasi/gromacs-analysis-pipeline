# GROMACS Post-Simulation Analysis Pipeline

This repository contains a unified and reusable Bash script to automate comprehensive post-simulation analysis of GROMACS molecular dynamics trajectories. The pipeline integrates trajectory preprocessing, essential structural analyses, and advanced techniques like Principal Component Analysis (PCA) and Free Energy Landscape (FEL) mapping.

Originally assembled by consolidating command histories across multiple simulation projects, this workflow streamlines repetitive tasks into a single robust script.

## 🔧 Key Features

- **Trajectory Pre-processing**: Removes periodic boundary conditions (PBC), centers the system, and prepares the trajectory for analysis.
- **RMSD**: Calculates backbone Root Mean Square Deviation across time to assess structural stability.
- **RMSF**: Calculates C-alpha Root Mean Square Fluctuation to identify flexible regions.
- **Radius of Gyration**: Assesses protein compactness over the simulation period.
- **Hydrogen Bond Analysis**: Counts and tracks intra-molecular hydrogen bonds.
- **Secondary Structure (DSSP)**: Tracks the evolution of secondary structural elements using DSSP.

### PCA (Principal Component Analysis)
- Constructs the covariance matrix of atomic displacements.
- Projects the trajectory along principal eigenvectors to capture dominant motions.

### FEL (Free Energy Landscape)
- Projects motion onto selected principal components.
- Generates a 2D Gibbs free energy surface to visualize metastable states.

- **Structure Extraction**: Exports snapshots (PDB) at defined time points for visualization or further analysis.

---

## ⚙️ Requirements

- **GROMACS** (tested with version 2022)
- **DSSP** (`dssp` binary should be in your system's `PATH`)
- **Perl**
- **Python 2.7**

**Helper Scripts** (must be in the same directory as the main script):
- `sham.pl`: Prepares data for `gmx sham`
- `xpm2txt.py`: Converts `.xpm` FEL maps to human-readable text

---

## 🚀 Usage

### Step 1: Prepare Your Directory
- Place `analysis_pipeline.sh` in a new directory.
- Copy your simulation files (`.tpr`, `.xtc`) into the same directory.
- Add `sham.pl` and `xpm2txt.py` alongside the main script.

### Step 2: Configure the Script
- Open `analysis_pipeline.sh` in a text editor.
- In the **USER CONFIGURATION** section, set the `PREFIX` variable to match your file base name.  
  For example, if your files are `system.tpr` and `system.xtc`, set `PREFIX="system"`.
- Modify other variables as needed (residue ranges, PC selections, etc.).

### Step 3: Make the Script Executable and Run

```bash
chmod +x analysis_pipeline.sh
./analysis_pipeline.sh

