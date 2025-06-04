#!/bin/bash
set -e

# Example script for performing unfolding

# Generate the supercell kpoint and easyunfold.json file
easyunfold generate -y --no-expand --out-file=no_expand.json Si/POSCAR Si_super_deformed/POSCAR Si/KPOINTS_band 

# Copy the supercell kpoints to the working directory
cp KPOINTS_no_expand Si_super_deformed_no_expand

# Run VASP to obtain the wave function at requested kpoints
cd Si_super_deformed_no_expand
bash run.sh
cd ../

# Calculate the spectral weights and generate effective band structure plots
easyunfold unfold --data-file=no_expand.json calculate Si_super_deformed_no_expand/WAVECAR
easyunfold unfold --data-file=no_expand.json  plot --emax=10 --emin=-10 --out-file=unfold_no_expand.png --intensity=2
