# MgO with atomic projections

:::{note}
The files needed for this example are provided in the 
[examples/MgO](https://github.com/SMTG-UCL/easyunfold/tree/main/examples/MgO) folder.
:::

Often it is useful to know the various contributions of atoms in the structure to the electronic bands 
in the band structure, to analyse the chemistry and orbital interactions at play in the system. This 
can be computed for unfolded bands as well.

For a normal band structure calculation, the contributions can be inferred by colouring the band 
according to the elemental contributions, which can be done using [sumo](https://github.com/SMTG-UCL/sumo).

```{figure} ../../examples/MgO/MgO/band.png
:width: 400px
:alt: Band structure of MgO primitive cell.

Band structure of MgO with atomic contribution </figcaption>
```

Similar plots can be generated for unfolded band structures. However, because the unfolded spectral 
function itself contains both the *location* of the band and its *intensity*, adding a third 
dimension of information (atomic projection) can be tricky to visualise.

In this example, we unfold the bands from a MgO 2x1x2 supercell with a Mg atom displaced to break 
symmetry. The procedure is essentially the same as described in the 
[Si supercell example](https://smtg-ucl.github.io/easyunfold/examples/example_si222.html).

The only difference here is that we turn on the calculation of orbital projections in `VASP` with 
`LORBIT = 11` (`12`, `13` and `14` will also work) in the `INCAR` file, and then use the `plot-projections` subcommand 
when plotting the unfolded band structure:

```bash
easyunfold unfold plot-projections --procar MgO_super/PROCAR --atoms="Mg,O" --combined --emin=-6 \
--emax=20 --intensity 6.5
```

Note that the path of the `PROCAR` is passed along with the desired atom projections (`Mg` and `O` here). 

:::{tip}
If the _k_-points have been split into multiple calculations (e.g. hybrid DFT band structures), the `--procar` option 
should be passed multiple times to specify the path to each split `PROCAR` file (i.e. 
`--procar calc1/PROCAR --procar cal2/PROCAR ...`).
:::

:::{note}
The atomic projections are not stored in the `easyunfold.json` data file, so the `PROCAR` file(s) should be kept for 
replotting in the future.
:::

The `--combined` option creates a combined plot with different colour maps for each atomic grouping.
The spectral intensity is used to define the transparency (`alpha`) allowing the fusion of multiple
projections into a single plot.

```{figure} ../../examples/MgO/unfold_project.png
:alt: MgO unfolded band structure with atomic projections
:width: 600px
:align: center
Unfolded MgO band structure with atomic projections. 
```

The default colour map for atom projections is red, green, blue and purple for the 1st, 2nd, 3rd and 4th atomic species 
specified. This can be changed with the `--colours` option, as well as several other plotting/visualisation settings – 
see the output of `easyunfold unfold plot-projections --help` for more details. If we wanted to plot the atomic 
projections with the same colouring scheme as the `sumo` plot above (i.e. red for Mg and blue for O), we can use:

```bash
easyunfold unfold plot-projections --procar MgO_super/PROCAR --atoms="Mg,O" --combined --emin=-6 \
--emax=20 --intensity 6.5 --colours "r,b"
```

```{figure} ../../examples/MgO/unfold_project_rb.png
:alt: MgO unfolded band structure with specified colours
:width: 600px
:align: center
Unfolded MgO band structure with atomic projections; red for Mg and blue for O atoms. 
```

:::{note}
In order to specify the atomic projections with `--atoms`, the `POSCAR` or `CONTCAR` file for the supercell must be 
present. If this is not the case, or we want to use projections from only specific atom subsets in the supercell, we can 
alternatively use the `--atoms-idx` tag. This takes a string of the form `a-b|c-d|e-f` where `a`, `b`, `c`, `d`, `e` and
`f` are integers corresponding to the atom indices in the VASP structure file (i.e. `POSCAR`/`CONTCAR`, corresponding 
to the `PROCAR` being used to obtain the projections). Different groups are separated by `|`, and `-` 
can be used to define the range for each projected atom type. A comma-separated list can also be used instead of ranges 
with hyphens. Note that 1-based indexing is used for atoms, matching the convention in VASP, which is then converted to 
zero-based indexing internally in python. In this example, we could set `--atoms-idx="1-4|5-8"` to get the same result
as `--atoms="Mg,O"` (but without the figure legend).
:::

In some cases, especially if there are many projection to be plotted at the same time, it can be clearer to create 
separated plots for each. This is the default behaviour for `plot-projections`, when `--combined` is not specified.

```bash
easyunfold unfold plot-projections --procar MgO_super/PROCAR --atoms="Mg,O" --emin=-6 \
--emax=20 --intensity 6.5
```

```{figure} ../../examples/MgO/unfold_project_sep.png
:width: 800 px
:alt: Projected MgO band structure  

Unfolded MgO band structure with atomic projections plotted separately.
```

:::{tip}
There are _many_ customisation options available for the plotting functions in `easyunfold`. See `easyunfold plot -h` or 
`easyunfold unfold plot-projections -h` for more details!
:::
