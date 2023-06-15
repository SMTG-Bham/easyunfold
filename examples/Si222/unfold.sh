#!/bin/bash
set -e

# Example script for performing unfolding

# Generate the supercell kpoint and easyunfold.json file
easyunfold generate Si/POSCAR Si_super_deformed/POSCAR Si/KPOINTS_band -y
# Copy the supercell kpoints to the working directory
cp KPOINTS_easyunfold Si_super_deformed

# Run VASP to obtain the wave function at requested kpoints
cd Si_super_deformed
bash run.sh
cd ../

# Calculate the spectral weights and generate effective band structure plots
easyunfold unfold calculate Si_super_deformed/WAVECAR
easyunfold unfold plot --emax=10 --emin=-10
