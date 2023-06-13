# Unfolding a Si 2x1x1 supercell using CASTEP

Below is a step-by-step guide for unfolding a 2x1x1 perfect supercell sampled with a simple $\Gamma - X$ path using CASTEP.

## Generate the project file and kpoints for supercell calculation

:::{note} 
The files needed as provided in expamples/Si222-castep. This guide assumes the current
working directory is located at the root of that folder.
:::

First, generate the supercell kpoints:

```bash
easyunfold generate Si_prim.cell Si_211.cell band.cell --code castep
```

:::{tip}
Here, `band.cell` contains the band structure path and the labels of the high symmetry points. The [sumo](https://github.com/SMTG-UCL/sumo) package can be used to generate it with ease.
:::

The command above also generates an  `easyunfold.json` file in the current direction containing information about the unfolding.
Name of this file can be modified with the `--out-file` command-line argument.

Information stored in this file can be inspected with command:

```bash
$ easyunfold unfold status
Loaded data from easyunfold.json
Supercell cell information:
        Space group number: 227
        Internation symbol: Fd-3m
        Point group: m-3m

Primitive cell information:
        Space group number: 227
        Internation symbol: Fd-3m
        Point group: m-3m

No. of k points in the primitive cell           : 21
No. of expanded kpoints to be calculated cell   : 60 (61)
No. of rotations in the primitive cell          : 48
No. of rotations in the super cell              : 24

Path in the primitive cell:
   \Gamma    : 1    
   X         : 21   
```

In this example, the path goes from $\Gamma$ to $X$.

Copy the generated `cell` file to the supercell calculation folder:

```bash
cp easyunfold_sc_kpoints.cell Si_211_unfold/
```


## Performing supercell band structure calculation

Band structure calculation in CASTEP is done by having the following keys in the `<seed>.param` file:

```
task : spectral
spectral_task: bandstructure
write_orbitals : true
```

The last key is necessary for band unfolding, as it enables the wave function of the band structure calculation to be written to the disk in the `<seed>.orbitals` file.

Assuming CASTEP is installed, the following command will launch the calculation:

```bash
cd Si_211_unfold
mpirun -np 4 castep.mpi easyunfold_sc_kpoints 
```

CASTEP will first perform a standard self-consistent field loop to find the ground state electron density (and hence the potential).
Afterwards, (in the same run) it switches to spectral calculation and solve for the eigenvalues at each requested k point with fixed potential (density).


## Perform unfolding

We now use calculate the spectral weights and also record the VBM:

```
cd ../
easyunfold unfold calculate Si_211_unfold/easyunfold_sc_kpoints.orbitals
```

Plot the unfolded band structure:

```bash
easyunfold unfold plot
```


Output:

```{figure} ../../examples/Si-castep/unfold.png
:alt: Spectral function
:width: 600px

Spectral function of the unfolded bands.
```


## Band structures in the primitive and the super cell

Below is the band structure of the primitive cell, along the same path:

```{figure} ../../examples/Si-castep/Si_band/band.png
:alt: Primitive cell band structure
:width: 600px

Primitive cell band structure of Si.
```

and the unfold band structure resembles it well.
This is not surprising, since a perfect supercell is used.

Below is the band structure in the supercell, blindly following the same path:

```{figure} ../../examples/Si-castep/Si_211_band/band.png
:alt: Primitive cell band structure
:width: 300px

Supercell cell band structure of Si in a 2x1x1 supercell
```

which, as expected, is the same as the primitive cell if it is folded back along the mid point between $\Gamma$ and $X$.