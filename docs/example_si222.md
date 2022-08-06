## Unfolding 2x2x2 Si supercell with a displaced atom

### Generate project file and kpoints for supercell calculation

!!! note 
    
    The files needed as provided in expamples/Si222. This guide assumes the current
    working directory is located at the root of that folder.

First, generate the supercell kpoints:

```
easyunfold generate Si/POSCAR Si_super_deformed/POSCAR Si/KPOINTS_band
```

Copy the kpoints to the supercell calculation folder:

```
cp KPOINTS_easyunfold Si_supercell_deformed
```

This generates an  `easyunfold.json` file in the current direction containing information about the unfolding.
Name of this file can be modified with the `--out-file` commandline argument.

Information stored in this file can be inspected with command:

```
$ easyunfold unfold status
Loaded data from easyunfold.json
Primitive cell information:
        Space group number: 160
        Internation symbol: R3m
        Point group: 3m

Supercell cell information:
        Space group number: 227
        Internation symbol: Fd-3m
        Point group: m-3m

No. of k points in the primitive cell           : 73
No. of expanded kpoints to be calculated cell   : 132 (136)
No. of rotations in the primitive cell          : 48
No. of rotations in the super cell              : 6

Path in the primitive cell:
   \Gamma    : 1    
   L         : 21   
   W         : 38   
   X         : 50   
   \Gamma    : 73   
Please run the supercell band structure calculation and run `unfold calculate`.
```

### Performing supercell band structure calculation

Band structure calculation in VASP normally involves two steps. First, a normal single point calculation is performed to obtain the charge density.
Afterwards, a none self-consistent calculation is carried out to compute the eigenvalues of the kpoints along the band structure paths defined.

First, ensure the kpoints for SCF is used and run the supercell calculation. The `ICHARG=11` must commented out in the INCAR for the single point calculation:

```
cd Si_supercell_deformed
cp KPOINTS_scf KPOINTS
sed -i 's/^ICHARG = 11/!ICHARG = 11/g' INCAR
mpirun -np 4 vasp_std 
```

Now run the band structure calculation with `ICHARG=11`, and the kpoints mapped to the supercell from the primitive cell path:

```
sed -i 's/.*ICHARG = 11/ICHARG = 11/g' INCAR
cp KPOINTS_easyunfold KPOINTS
mpirun -np 4 vasp_std 
```

### Perform unfolding

Calculate the weights and record the VBM:

```
cd ../
easyunfold unfold calculate Si_super_deformed/WAVECAR
```

!!! note 

    If you don't wnat to run the VASP calculation by yourself, the calculated `WAVECAR` and `vasprun.xml` for this example with:

    ```
    wget -o Si_super_deformed/WAVECAR https://www.dropbox.com/s/3cmn2epw7d290jd/WAVECAR?dl=1
    wget -o Si_super_deformed/vasprun.xml https://www.dropbox.com/s/ota78qqdvxsehmi/vasprun.xml?dl=1
    ```

Plot the unfolded band structure:

```
easyunfold unfold plot
```


Output:

<figure markdown>
  ![Spectral function](img/Si_222_unfold.png){ width="400" }
  <figcaption> Spectral function of the unfolded bands </figcaption>
</figure>


Band structure of the primitive cell:

<figure markdown>
  ![Primitive cell band structure](img/Si_222_band.png){ width="400" }
  <figcaption> Primitive cell band structure of Si </figcaption>
</figure>


Note the appearance of extra branches due to symmetry breaking.  

### What happens if the symmetry is not taken into account?


We can create a new unfolding project (data file) using the following command:

```
easyunfold generate Si/POSCAR Si_super_deformed/POSCAR Si/KPOINTS_band --no-expand --out-file no-expand.json
```

Swap the `KPOINTS` to the new file

```
cp KPOINTS_no-expand Si_super_deformed/KPOINTS
cd Si_super_deformed
mpirun -np 4 vasp_std
cd ../
easyunfold --help no-expand.json calculate Si_super_deformed/WAVECAR
easyunfold unfold --data-file  no-expand.json  plot --out-file unfold_no-expand.png
```

output:

<figure markdown>
  ![Spectral function](img/Si_222_unfold_no-expand.png){ width="400" }
  <figcaption> Spectral function of the unfolded bands </figcaption>
</figure>

Comparing with the one above, there are breaking of the bands and some branches are missing (near the gamma point).


Nevertheless, by not expanding the kpoint paths, fewer supercell kpoints need to be calculated. 
 
```
$ easyunfold unfold --data-file  no-expand.json  plot --out-file unfold_no-expand.png

Loaded data from no-expand.json
Using a reference energy of 5.284 eV
Unfolded band structure saved to unfold_no-expand.png

$ easyunfold unfold --data-file  no-expand.json  status

Loaded data from no-expand.json
Primitive cell information:
        Space group number: 160
        Internation symbol: R3m
        Point group: 3m

Supercell cell information:
        Space group number: 227
        Internation symbol: Fd-3m
        Point group: m-3m

No. of k points in the primitive cell           : 73
No. of expanded kpoints to be calculated cell   : 70 (73)
No. of rotations in the primitive cell          : 48
No. of rotations in the super cell              : 6

Path in the primitive cell:
   \Gamma    : 1    
   L         : 21   
   W         : 38   
   X         : 50   
   \Gamma    : 73   
Unfolding had been performed - use `unfold plot` to plot the spectral function.
```