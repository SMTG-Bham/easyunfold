cp KPOINTS_scf KPOINTS
rm WAVECAR CHGCAR CHG
sed -i 's/.*ICHARG.*/ICHARG = 2/g' INCAR
echo "Running SCF"
mpirun -np 4 vasp_std
sed -i 's/.*ICHARG.*/ICHARG = 11/g' INCAR
cp KPOINTS_no_expand KPOINTS
echo "Running BS"
mpirun -np 4 vasp_std
