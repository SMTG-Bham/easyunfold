#!/bin/bash
set -e

# Example script for performing unfolding

# Generate the supercell kpoint and easyunfold.json file
easyunfold generate MgO/POSCAR MgO_super/POSCAR MgO/KPOINTS_band -y
# Copy the supercell kpoints to the working directory
cp KPOINTS_easyunfold MgO_super

# Run VASP to obtain the wave function at requested kpoints
cd MgO_super
bash run.sh
cd ../

# Calculate the spectral weights and generate effective band structure plots
easyunfold unfold calculate MgO_super/WAVECAR
easyunfold unfold plot --emax=20 --emin=-15
easyunfold unfold plot-projections --emax=20 --emin=-15 --atoms-idx="1-4|5-8" --out-file unfold_project.png --procar MgO_super/PROCAR
