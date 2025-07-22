# GROMACS Post-Simulation Analysis Pipeline

This repository contains a unified and reusable Bash script to automate comprehensive post-simulation analysis of GROMACS molecular dynamics data. The pipeline handles tasks like trajectory cleaning, RMSD, RMSF, secondary structure analysis (DSSP), and advanced analyses such as Principal Component Analysis (PCA) and Free Energy Landscape (FEL) generation.

This script was collated from a series of command history files from multiple simulation projects to create a single, robust workflow.

## Key Features

The pipeline automates the following analyses, which were common across the provided project files:
- [cite_start]**Trajectory Pre-processing**: Removes periodic boundary conditions (PBC) and centers the protein in the box[cite: 559, 662, 865].
- [cite_start]**Root Mean Square Deviation (RMSD)**: Calculates the RMSD of the protein backbone over time[cite: 33, 352, 359, 511, 566, 669, 779, 872, 879].
- [cite_start]**Root Mean Square Fluctuation (RMSF)**: Calculates the RMSF of C-alpha atoms to identify flexible regions[cite: 40, 573, 677, 879].
- [cite_start]**Radius of Gyration**: Measures the compactness of the protein over the simulation[cite: 124, 365, 550, 635, 769, 994].
- [cite_start]**Hydrogen Bond Analysis**: Calculates intra-protein hydrogen bonds[cite: 340, 473, 505].
- [cite_start]**Secondary Structure Analysis**: Uses DSSP to analyze the evolution of secondary structure elements[cite: 192, 369, 553, 639, 840, 1058].
- **Principal Component Analysis (PCA)**:
    - [cite_start]Calculates the covariance matrix[cite: 55, 199, 375, 582, 686, 888].
    - [cite_start]Analyzes eigenvectors and projects the trajectory onto the principal components[cite: 68, 78, 90, 210, 223, 235, 247, 261, 389, 402, 413, 424, 592, 603, 616, 696, 708, 721, 735, 898, 907, 920, 933].
- **Free Energy Landscape (FEL)**:
    - [cite_start]Generates a 2D projection of the trajectory onto two principal components[cite: 103, 301, 455, 747, 945].
    - [cite_start]Creates a Gibbs free energy landscape from the projection[cite: 112, 309, 464, 755, 955].
- [cite_start]**Structure Extraction**: Dumps PDB files of the protein at specific time points for visualization[cite: 44, 48, 143, 147, 152, 158, 162, 168, 172, 178, 182, 187, 274, 279, 284, 288, 292, 296, 311, 319, 436, 441, 446, 450, 576, 681, 760, 764, 786, 790, 796, 800, 806, 810, 815, 820, 824, 829, 835, 882, 958, 963, 967, 973, 984, 988, 1010, 1016, 1020, 1025, 1030, 1034, 1039, 1043, 1048, 1053].

## Requirements
- **GROMACS** (workflow based on version 2022)
- **DSSP**: The `dssp` executable must be in your system's PATH.
- **Perl**
- **Python 2.7**
- **Helper Scripts**: The following scripts (used in the original analysis) should be in the same directory as the pipeline script:
    - `sham.pl` (to prepare data for `gmx sham`)
    - `xpm2txt.py` (to convert FEL xpm files to text)

## Usage

1.  **Prepare Your Directory**:
    - Place the `analysis_pipeline.sh` script in a new directory.
    - Copy your GROMACS input files (`.tpr` and `.xtc`) into the same directory.
    - Add the required helper scripts (`sham.pl`, `xpm2txt.py`).

2.  **Configure the Script**:
    - Open `analysis_pipeline.sh` in a text editor.
    - In the **USER CONFIGURATION** section, change the `PREFIX` variable to match the base name of your input files. For example, if you have `my_system.tpr` and `my_system.xtc`, set `PREFIX="my_system"`.
    - Adjust other variables like residue ranges or selected PCs for FEL as needed.

3.  **Make the Script Executable**:
    ```sh
    chmod +x analysis_pipeline.sh
    ```

4.  **Run the Pipeline**:
    ```sh
    ./analysis_pipeline.sh
    ```
    The script will execute all analysis steps sequentially, creating output files prefixed with your chosen `PREFIX`.

## Workflow Overview

The script is organized into logical sections:

1.  **Trajectory Pre-processing**: Creates a clean, whole trajectory (`_noPBC.xtc`) for analysis.
2.  **Standard Analysis**: Calculates RMSD, RMSF, Radius of Gyration, and Hydrogen Bonds.
3.  **Secondary Structure**: Runs DSSP and generates a plot of secondary structure vs. time.
4.  **Principal Component Analysis (PCA)**: Calculates and analyzes the principal modes of motion.
5.  **Free Energy Landscape (FEL)**: Generates a 2D FEL plot based on the most significant principal components.
6.  [cite_start]**Optional Analyses**: Contains commented-out code blocks for specialized analyses, such as domain-specific RMSD [cite: 518, 526, 534, 542] [cite_start]or analyzing water molecules near a binding site[cite: 478, 484, 490, 492, 495, 499, 502, 505]. These can be enabled by uncommenting the code and configuring the residue selections.
7.  **Structure Extraction**: Saves PDB snapshots at regular intervals for visualization.

## Simulation Setup Commands

The provided history files also contained commands for running the simulations themselves. While not part of the analysis pipeline, a typical CHARMM-GUI-based workflow includes these steps:

- [cite_start]**Minimization**: `gmx grompp ... && gmx mdrun ...`[cite: 1, 644, 845].
- [cite_start]**Equilibration**: `gmx grompp ... && gmx mdrun ...`[cite: 9, 652, 854].
- [cite_start]**Production**: `gmx grompp ... && gmx mdrun ...`[cite: 17, 22, 659, 862].

These are generally run once to generate the trajectory, which is then analyzed by the main pipeline script.
