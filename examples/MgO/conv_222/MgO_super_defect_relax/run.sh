cp KPOINTS_scf KPOINTS
rm WAVECAR CHGCAR CHG
exe="mpirun -np 4 singularity run ~/appdir/vasp-6.3.0-fftw-openblas.sif  vasp_std"
sed -i 's/.*ICHARG.*/ICHARG = 2/g' INCAR
echo "Running SCF"
eval $exe 
sed -i 's/.*ICHARG.*/ICHARG = 11/g' INCAR
cp KPOINTS_easyunfold KPOINTS
echo "Running BS"
eval $exe 
