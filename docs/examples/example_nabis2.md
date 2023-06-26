# Cation-Disordered NaBiS<sub>2</sub> with atomic projections

:::{note}
The files needed for reproducing this example are provided in the 
[examples/NaBiS2](https://github.com/SMTG-UCL/easyunfold/tree/main/examples/NaBiS2) folder. 
Note that the `PROCAR.gz` file will need to be decompressed with `gzip -d PROCAR.gz` if recalculating 
and reproducing these example plots.
:::

In this example, we unfold the bands from a 80-atom special-quasirandom (SQS) supercell of 
NaBiS<sub>2</sub>, where the Na and Bi cations are quasi-randomly distributed, in order to simulate 
randomised cation disorder in the material. 

:::{tip}
SQS supercells can be generated using tools like 
[icet](https://icet.materialsmodeling.org/moduleref_icet/tools.html#module-icet.tools.structure_generation)
or [ATAT](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).
:::

These results were published in Y. T. Huang & S. R. Kavanagh et al. 2022 [^1], and an early version of 
`easyunfold` was used for the similar AgBiS$_2$ in Y. Wang & S. R. Kavanagh et al. 2022 [^2], 
with these plots demonstrating the key 
differences in electronic structure and thus photovoltaic performance between these two materials.

We have previously calculated the `easyunfold.json` file from the calculation using `easyunfold calculate WAVECAR`. 
Using the default plotting options with `easyunfold unfold plot`, we obtain the following unfolded band structure:

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot.png
:alt: NaBiS2 unfolded band structure
:width: 400px
Unfolded band structure of $\ce{NaBiS2}$
```

This is nice, but we can make the plot a little clearer by adjusting some of the parameters like the intensity scaling
(via `--vscale`) and the colour map (via `--cmap`). Below we set `--vscale 0.4` to increase the colourmap intensity, 
and try `BuPu`, `viridis` and `bone_r` from left to right below: 

```bash
easyunfold unfold plot --vscale 0.4 --cmap "BuPu"
easyunfold unfold plot --vscale 0.4 --cmap "viridis"
easyunfold unfold plot --vscale 0.4 --cmap "bone_r"
```

BuPu             |  viridis        |  bone_r
:-------------------------:|:-------------------------:|:-------------------------:
![](../../examples/NaBiS2/NaBiS2_unfold-plot_BuPu.png)  |  ![](../../examples/NaBiS2/NaBiS2_unfold-plot_viridis.png) |  ![](../../examples/NaBiS2/NaBiS2_unfold-plot_bone_r.png)

We can also plot the unfolded band structure with atomic projections, which is useful for understanding the electronic 
structure of the material. In this case, we are curious as to which atoms are contributing to the band edges, and so 
the atomic projections will be useful. For this, we need the `PROCAR` output from VASP with the atomic and orbital 
projection information, and so `LORBIT` should be set to `11` or `12` in the `INCAR` for the original calculation.

When plotting the unfolded band, the `plot-projections` subcommand is used with the `--procar <PROCAR>` and 
`--atoms-idx <atoms-idx>` options:

```bash
easyunfold unfold plot-projections --atoms-idx="1-20|21-40|41-80" --procar PROCAR  --intensity=2  --combined
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj.png
:alt: NaBiS2 atom-projected unfolded band structure
:width: 400px

Unfolded band structure of NaBiS<sub>2</sub> with atomic contributions.
```

From this plot, we can see that Sulfur anions (in blue) dominate the valence band, while Bismuth cations (in green) 
dominate the conduction band, with minimal contributions from the Sodium cations as expected.

Note that the path of the `PROCAR` is passed along with the group of atoms.
In this example, the first 20 atoms are `Na`, the second 20 are `Bi` and the last 40 are `S`. Different groups are 
separated by `|`, and `-` can be used to define the range.
Note that we use "1-based indexing" for the atoms here, matching the VASP format (i.e. the index of the first atom is 1, 
not zero as in python).

:::{note}
The atomic projections are not stored in the `easyunfold.json` data file, so the `PROCAR` file should be 
kept for replotting in the future.
:::

While the main conclusions of S dominating the valence band and Bi dominating the conduction band are clear from the 
plot above, the high intensity of the S projections in the valence band makes the Bi contributions in the conduction 
band more faint and less apparent. 

So, we can create separated plots for each atom type to make their individual contributions more clear:

```bash
easyunfold unfold plot-projections --atoms-idx="1-20|21-40|41-80" --procar PROCAR  --cmap="bone_r" --vscale 0.4
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_sep.png
:alt: NaBiS2 atom-projected separated unfolded band structure 
:width: 800px

Unfolded band structure of NaBiS<sub>2</sub> with atomic contributions plotted separately.
```


An alternative option here is also to just plot only the contributions of Na (`1-20`) and Bi (`21-40`) 
cations, with no S projections:
```bash
easyunfold unfold plot-projections --atoms-idx="1-20|21-40" --procar PROCAR  --intensity=2  --combined --colors="r,g"
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_noS.png
:alt: NaBiS2 atom-projected unfolded band structure, no S 
:width: 400px

Unfolded band structure of NaBiS<sub>2</sub> with atomic contributions of only Na and Bi.
```

While this plot isn't the most aesthetic, it clearly shows that Bi (green) contributes to both the conduction band and 
(less so) valence band states, but Na (red) doesn't contribute significantly at near the band edges 
(it's a spectator ion!). 

[^1]: [Huang, YT., Kavanagh, S.R., Righetto, M. et al. Strong absorption and ultrafast localisation in NaBiS2 nanocrystals with slow charge-carrier recombination. Nat Commun 13, 4960 (2022)](https://www.nature.com/articles/s41467-022-32669-3) 
[^2]: [Wang, Y., Kavanagh, S.R., Burgués-Ceballos, I. et al. Cation disorder engineering yields AgBiS2 nanocrystals with enhanced optical absorption for efficient ultrathin solar cells. Nat. Photon. 16, 235–241 (2022).](https://www.nature.com/articles/s41566-021-00950-4)