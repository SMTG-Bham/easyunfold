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
`LORBIT = 11` in the `INCAR` file, and then use the `plot-projections` subcommand when plotting the 
unfolded band structure:

```bash
easyunfold unfold plot-projections --procar MgO_super/PROCAR --atoms-idx="1-4|5-8" \ 
--out-file unfold_project.png --combined --emin=-15 --emax=15
```

Note that the path of the `PROCAR` is passed along with the desired groupings of atoms (with 
`--atoms-idx`).
In this example, the first four atoms are `Mg` (`1-4`) and the last four are `O` (`5-8`), and we would
like to show the elemental contributions to the bands. Different groups are separated by `|`, and `-` 
can be used to define the range. Note that use of one-based indexing for atoms, although in python 
zero-based indexing is used internally.

:::{note}
The atomic projections are not stored in the `easyunfold.json` data file, so the `PROCAR` file should be 
kept for replotting in the future.
:::

The `--combined` option creates a combined plot with different colour maps for each atomic grouping.
The spectral intensity is used to define the transparency (`alpha`) allowing the fusion of multiple
projections into a single plot.

```{figure} ../../examples/MgO/unfold_project.png
:alt: MgO unfolded band structure
:width: 400px
:align: center
Unfolded MgO band structure with atomic projections. Red for Mg and green for O atoms. 
```


In some cases, especially if there are many projection to be plotted at the same time, it can be clearer to create separated plots for each.

```bash
easyunfold unfold plot-projections  --procar MgO_super/PROCAR --atoms-idx="1-4|5-8" \
--out-file unfold_project_sep.png --emax=22 --emin=-18
```

```{figure} ../../examples/MgO/unfold_project_sep.png
:width: 800 px
:alt: Projected MgO band structure  

Unfolded MgO band structure with atomic projections plotted separately (Mg left and O right).
```
