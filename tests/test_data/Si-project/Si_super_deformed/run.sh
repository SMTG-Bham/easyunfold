cp KPOINTS_scf KPOINTS
rm WAVECAR CHGCAR CHG
sed -i 's/.*ICHARG.*/ICHARG = 2/g' INCAR
echo "Running SCF"
vasp_std
sed -i 's/.*ICHARG.*/ICHARG = 11/g' INCAR
cp KPOINTS_easyunfold KPOINTS
echo "Running BS"
vasp_std
