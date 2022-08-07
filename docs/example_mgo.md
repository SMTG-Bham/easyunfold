## Unfolding MgO band structure with atomic projections

!!! note
    
    Relevant files can be found in the `examples/MgO` folder.

In some cases, it is useful to know the atomic contributions of the bands.
This can be done for unfolded bands as well.

For a normal band structure, the contributions can be inferred by colouring the band according to the elemental contributions.

<figure markdown>
  ![MgO band structure](img/MgO_band.png){ width="400" }
  <figcaption> Band structure of MgO with atomic contribution </figcaption>
</figure>


Similar plots can be generate for unfolded band structure. However, because the spectral function itself contains both the *location* of the band and its *intensity*, adding a third information regarding the projection can be tricky.

In this example, we unfold the bands from a MgO 2x1x2 supercell with the first Mg atom displaced. The procedure is essentially the same as the Si supercell example.

The only difference is that we turn on the calculation of orbital projections with `LORBIT=11` in the `INCAR` file.

When plotting the unfolded band, the `plot-projections` subcommand is used:

```
easyunfold unfold plot-projections  --procar MgO_super/PROCAR \ 
--atoms-idx="1-4|5-8" --out-file unfold_project.png --combined --cmap="Greens|Reds"
```

Note that the path of the `PROCAR` is passed along with the group of atoms.
In this example, the first four atoms are `Mg` the last four are `O`, and we would
like to show the contribution of the band based on the elements.
Different groups are separated by `|`, and `-` can be used to define the range.
Note that use of one-based indexing for atoms, although in python zero-based indexing is used internally.


!!! note

    The projections are not stored in the `easyunfold.json` data file. So the `PROCAR` is should be kept for replotting in the future.

The `--combined` option creates a combined plot with each group colour with different colour maps.
The spectral intensity is used to define the transparency (`alpha`) allowing the fusion of multiple
projections into a single plot.

<figure markdown>
  ![MgO unfolded band structure](img/MgO_unfold_project.png){ width="400" }
  <figcaption> Unfolded MgO band structure with projections. Green for Mg atoms and red for O atoms. </figcaption>
</figure>

In some cases, especially if there are many projection to be plotted at the same time, it can be clearer to create separateted plots for each.

```
easyunfold unfold plot-projections  --procar MgO_super/PROCAR --atoms-idx="1-4|5-8" --out-file unfold_project_sep.png
```

<figure markdown>
  ![MgO unfolded band structure](img/MgO_unfold_project_sep.png){ width="800" }
  <figcaption> Unfolded MgO band structure with projections plotted separately. </figcaption>
</figure>






