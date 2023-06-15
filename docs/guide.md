# Guide

Main goal of this tool is to make the unfolding process easier.
To generate a unfolded band structure, one typically needs to perform the following step:

1. Create a primitive cell, and generate a k point path for this primitive cell.
2. Create a supercell, and obtain its optimised structure.
3. Generate a series of kpoints in the supercell to be calculated.
4. Perform a band structure calculation using the supercell, and save its wave function.
5. Run post-processing to obtain the unfolded band structure.

The supercell usually contains certain defects, or a special quasi random structure.
In both cases, its symmetry is lowered when compared to the perfect primitive cell.
Hence, for a given kpoint path in the primitive cell, additional kpoints may need to be sampled, and the extracted spectral weights need to be averaged in the end to obtained the effective band structure (EBS).


At the moment, two plane wave DFT codes, VASP[^vasp] and CASTEP[^castep], are supported. In principle, other PW code can be supported easily if the wavefunction can be read in.

Below is a guide for performing unfolding when using VASP.

## Step 1 - Generate the kpoints path of the primitive cell

This can be done by well established packages such as [seekpath](https://www.materialscloud.org/work/tools/seekpath).
Be careful that the "standardised" primitive cell may be different from input structure,
and the generated path is correct for the standard primitive cell only.
We recommand using [sumo](https://github.com/SMTG-UCL/sumo) for generating the kpoints, which provides a nice command line interface:

```
sumo-kgen -p POSCAR
```



Care should be taken if one uses the initial structure for further supercell generation, do verify that the lattice parameters are identical between the two.
A `POSCAR_prim` file will be written out if `sumo` think the primitive cell is different from the input structure.
The kpoints along the path is written to `KPOINTS_band`.

:::{tip}
`sumo` can be installed with `pip`:

```bash
pip install sumo
```

:::

## Step 2 - Generate the kpoints to be used for the supercell calculation

At this point, you should have your relaxed supercell structure (which may have a lower symmetry).
The set of kpoints for the supercell band structure can be generated with:

```
easyunfold generate primitive/POSCAR supercell/POSCAR primitive/KPOINTS_band --matrix "2 2 2"
```

for hybrid functional calculations, it an be useful to split the kpoints into multiple calculations for reduced costs or memory consumptions.

Note that the `--matrix` input is for setting the transformation matrix such that

```
cell_super = M @ cell_primitive
```

where `cell_super` and `cell_primitive` are (3,3) matrix made of **row vectors**.
If `M` is non-diagonal, all nine elements must be passed in a row-major order.

It is possible to omit `--matrix`` if the supercell is perfectly commensurate with the primitive cell.
This can be the case if the supercell calculation did not undergo cell relaxation.

If cell relaxation did take place, it is important to note that the unfolded band structure is not for exact  original primitive cell, but for a primitive cell deformed in a similar way as the supercell.

A `easyunfold.json` file will be written which contains the information of the unfolding.
The kpoints needed to be calculated for the supercell is written to a file named `KPOINTS_easyunfold`.
It is possible to change the name `easyunfold` by passing a explicit tag with the command line `--out-file`.

For hybrid functional calculations, you may want to split the kpoints into multiple calculations:

```
easyunfold generate primitive/POSCAR supercell/POSCAR primitive/KPOINTS_band --matrix "2 2 2" --scf-kpoints IBZKPT --nk-per-split 60
```

This will generate files named as `KPOINTS_easyunfold_001`, `KPOINTS_easyunfold_002`, each containing 60 kpoints.
If a `IBZKPT` file is the provided, its kpoints will be included with their original weights,  and all of the kpoints included by easyunfold will be zero-weighted.
This is necessary for hybrid functional calculations where the electronic minimisation must be conducted self-consistently (e.g. `ICHARG<10`).

:::{tip}

Input files for CASTEP can be generated with option `--code castep`. 
In this case, the `<seed>.cell` file should be passed containing the primitive cell kpoints that are
stored under the `spectral_kpoints_list` block.
The `cell` file of the supercell structure will be used as a template for generating a single `cell` file containing the 
supercell kpoints stored under the `spectral_kpoints_list` block.
The choice of the DFT code is stored in the data file (`easyunfold.json`) and will be used automatically the later steps.

:::


## Step 3 - Perform the unfolding

At this point, a supercell calculation should be completed with a `WAVECAR` (containing the plane wave coefficients) written containing all of the kpoints in the `KPOINTS_easyunfold` file generated.
This is typically a non self-consistent calculation with `ICHARG=11` for standard DFT, or a self-consistent calculation with zero-weighted kpoints if hybrid functional is used.

To compute the spectral weights, run the following command:

```
easyunfold unfold calculate WAVECAR
```

This command compute the spectral weight and save them into the  `easyunfold.json` file.
You can load the `easyunfold.json` file to read the spectral weights manually, or proceed with the command line interface to generate a plot.

:::{tip}

Wave functions (plane wave coefficients) in a CASTEP spectral calculation is stored in the `<seed>.orbitals` file,
and the `<seed>.check` file only contains the wave function for the self-consistent field calculation.
Note that the former is not written by default, and needs to be turned on by setting:

```
write_orbitals : true
```

in the `<seed>.param` file.

:::

If the kpoints has been split into multiple calculations (for example, for those involving hybrid functional), all of the wave function (e.g. `WAVECAR` for VASP) files need to be passed:

```
easyunfold unfold calculate calc1/WAVECAR calc2/WAVECAR
```

For large `WAVECAR`, it may take some time to parse and compute the weights.

## Step 4 - Plot the results

Simply do:

```
easyunfold unfold plot
```

to generate a plot of the spectral function.
It is possible to further customise the plot though command line arguments - see the help with `easyunfold unfold plot --help`.


[^vasp]: https://www.vasp.at
[^castep]: http://www.castep.org