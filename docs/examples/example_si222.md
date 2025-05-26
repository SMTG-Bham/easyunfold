# Silicon with a displaced atom

Below is a step-by-step guide for unfolding the electronic structure of a `2x2x2` supercell of 
crystalline silicon (Si) which contains a displaced atom, breaking symmetry.

## Generate the project file (`easyunfold.json`) and _k_-points for supercell calculation

:::{note} 
The files needed for this example are provided in the 
[examples/Si222](https://github.com/SMTG-Bham/easyunfold/tree/main/examples/Si222) folder. This 
guide assumes the current working directory is located at the root of that folder.
:::

First, generate the supercell _k_-points:

```bash
easyunfold generate Si/POSCAR Si_super_deformed/POSCAR Si/KPOINTS_band
```

Here, `KPOINTS_band` is the `KPOINTS` file corresponding to the band structure path for the primitive 
unit cell, which in this case was generated using `sumo-kgen` (see [Step 1](
https://smtg-Bham.github.io/easyunfold/guide.html#step-1-generate-the-kpoints-path-of-the-primitive-cell) 
of the tutorial docs page).

:::{note} 
In this example we are using GGA DFT, but if we were using hybrid DFT with `VASP`, we would need to use the `--scf-kpoints` option with `easyunfold generate`, as described in the tutorial docs page.
:::

This generates an  `easyunfold.json` file in the current direction containing information about the 
unfolding. The output name of this file can be modified with the `--out-file` commandline argument.

Information stored in this file can be inspected with the `easyunfold unfold status` command:

```bash
$ easyunfold unfold status
Loaded data from easyunfold.json
Primitive cell information:
        Space group number: 160
        International symbol: R3m
        Point group: 3m

Supercell cell information:
        Space group number: 227
        International symbol: Fd-3m
        Point group: m-3m

No. of k points in the primitive cell           : 73
No. of supercell k points                       : 103
No. of primitive cell symmetry operations       : 48
No. of supercell symmetry operations            : 6

Path in the primitive cell:
   \Gamma    : 1    
   L         : 21   
   W         : 38   
   X         : 50   
   \Gamma    : 73   
Please run the supercell band structure calculation and run `unfold calculate`.
```

Copy the _k_-points to the supercell calculation folder:

```bash
cp KPOINTS_easyunfold Si_supercell_deformed
```

## Perform the supercell band structure calculation

Semi-local (GGA) DFT band structure calculations in VASP normally involve two steps. First, a normal 
single point calculation is performed to obtain the self-consistent (SCF) charge density. Following 
this, a non-self-consistent calculation is carried out to compute the eigenvalues of the _k_-points 
along the defined band structure path.

First, we run our single point SCF supercell calculation, ensuring to use the appropriate converged SCF 
_k_-point mesh for the supercell. The `ICHARG = 11` tag must not be set or be commented out in the `INCAR` 
file for this single point calculation:

```bash
cd Si_supercell_deformed
cp KPOINTS_scf KPOINTS  # SCF kpoint mesh
sed -i 's/^ICHARG = 11/!ICHARG = 11/g' INCAR  # comment out ICHARG = 11
mpirun -np 4 vasp_std  # run the calculation
```

Now, with our converged SCF charge density, we run the GGA band structure calculation with `ICHARG = 11`, 
and the _k_-points mapped to the supercell from the primitive cell path:

```bash
sed -i 's/.*ICHARG = 11/ICHARG = 11/g' INCAR  # set ICARG = 11
cp KPOINTS_easyunfold KPOINTS  # supercell band structure kpoint path
mpirun -np 4 vasp_std  # run the calculation
```

Alternatively, there is a `unfold.sh` script in the 
[examples/Si222/Si_super_deformed](https://github.com/SMTG-Bham/easyunfold/tree/main/examples/Si222/Si_super_deformed) 
folder that can be used to perform these two steps above.

## Perform band unfolding

Unfold the supercell wave function (`WAVECAR`) and calculate the spectral weights:

```
cd ../
easyunfold unfold calculate Si_super_deformed/WAVECAR
```

:::{note} 
If you don't want to run the VASP calculation by yourself, the calculated `WAVECAR` and `vasprun.xml` 
for this example can be downloaded using git-lfs:


If git-lfs was not installed when you cloned the repository, install it via:

```
sudo apt install git-lfs
git lfs install
```

Download the files:

```
git lfs pull -I git lfs pull -I examples/Si222/Si_super_deformed/WAVECAR,examples/Si222/Si_super_deformed/vasprun.xml --exclude=""
```
:::

Plot the unfolded band structure:

```bash
easyunfold unfold plot
```


Output:

```{figure} ../../examples/Si222/unfold.png
:alt: Spectral function
:width: 400px

Spectral function of the unfolded bands.
```

:::{tip} 
See the [NaBiS<sub>2</sub> example](https://smtg-Bham.github.io/easyunfold/examples/example_nabis2.html) for tips on 
customising and prettifying the unfolded band structure plot. For example, you can use the `--intensity 3.5` option to increase the spectral function intensity.
:::

Note the appearance of extra branches compared to the band structure of the primitive cell (below), due 
to symmetry breaking from the displaced atom.

```{figure} ../../examples/Si222/band.png
:alt: Primitive cell band structure
:width: 400px

Primitive cell band structure of Si.
```

## What happens if symmetry is not properly taken into account?

It is quite common that the supercell has lower symmetry compared to the primitive cell. 
By default, `easyunfold` takes account of such symmetry breaking effect by including
additional _k_-points that no longer equivalent under the symmetry of the supercell cell.

In this example, we show what happens if we **do not** include the additional kpoints.
We can create a new unfolding project (`json` data file) using the following command:

```bash
easyunfold generate Si/POSCAR Si_super_deformed/POSCAR Si/KPOINTS_band --no-expand --out-file no-expand.json
```

Swap the `KPOINTS` to the new file, non-expanded `KPOINTS` file:

```bash
cp KPOINTS_no-expand Si_super_deformed_no_expand/KPOINTS
cd Si_super_deformed_no_expand
mpirun -np 4 vasp_std
cd ../
easyunfold unfold --data-file  no-expand.json calculate Si_super_deformed_no_expand/WAVECAR
easyunfold unfold --data-file  no-expand.json  plot --out-file unfold_no-expand.png
```

output:

```{figure} ../../examples/Si222/unfold_no-expand.png
:alt: Spectral function
:width: 400px

Spectral function of the unfolded bands without additional kpoints due to reduced symmetry.
```

Comparing this plot with the one above, we see that we get spurious band breaking (e.g. along $\Gamma - L$)
and some branches are missing (near the $\Gamma$ point).

Nevertheless, by not expanding the _k_-point paths, fewer supercell _k_-points need to be calculated:
 
```bash
$ easyunfold unfold --data-file  no-expand.json  plot --out-file unfold_no-expand.png

Loaded data from no-expand.json
Using a reference energy of 5.284 eV
Unfolded band structure saved to unfold_no-expand.png

$ easyunfold unfold --data-file  no-expand.json  status

Loaded data from no-expand.json
Supercell cell information:
        Space group number: 160
        International symbol: R3m
        Point group: 3m

Primitive cell information:
        Space group number: 227
        International symbol: Fd-3m
        Point group: m-3m

No. of k points in the primitive cell           : 73
No. of supercell k points                       : 70
No. of primitive cell symmetry operations       : 48
No. of supercell symmetry operations            : 6

Path in the primitive cell:
   \Gamma    : 1    
   L         : 21   
   W         : 38   
   X         : 50   
   \Gamma    : 73   
Unfolding has been performed - use `unfold plot` to plot the spectral function.
```

Note that in most cases one would always want to include the additional kpoints to correctly capture the effect of symmetry breaking.
The `--no-expand` option should be used with care and **only when there is no alternative**,
for example, 
when the expansion gives too many kpoints for very large supercells of special quasi-random structures.

:::{note}
One can always split the workload into multiple calculations with `--nk-per-split` to fit the computational resources available for individual calculations.
:::
