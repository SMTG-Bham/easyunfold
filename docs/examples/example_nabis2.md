# Disordered NaBiS<sub>2</sub> with atomic/orbital projections & DOS plotting

:::{note}
The files needed for reproducing this example are provided in the 
[examples/NaBiS2](https://github.com/SMTG-Bham/easyunfold/tree/main/examples/NaBiS2) folder.
:::

In this example, we unfold the bands from a 80-atom special-quasirandom (SQS) supercell of 
NaBiS<sub>2</sub>, where the Na and Bi cations are quasi-randomly distributed, in order to simulate 
randomised cation disorder in the material. 

:::{tip}
SQS supercells can be generated using tools like 
[icet](https://icet.materialsmodeling.org/moduleref_icet/tools.html#module-icet.tools.structure_generation)
or [ATAT](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).
:::

These results were published in Y. T. Huang & S. R. Kavanagh et al. 2022 [^1], and an early version of `easyunfold` was 
used for the similar AgBiS$_2$ in Y. Wang & S. R. Kavanagh et al. 2022 [^2], with these plots demonstrating the key 
differences in electronic structure and thus photovoltaic performance between these two materials.

## Standard Unfolded Band Structure
We have previously calculated the `easyunfold.json` file from the calculation using `easyunfold calculate WAVECAR`. 
Using `easyunfold unfold plot --intensity 2`, we obtain the following unfolded band structure:

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot.png
:alt: NaBiS2 unfolded band structure
:width: 400px
Unfolded band structure of $\ce{NaBiS2}$
```

### Visualisation Customisation: Colour map and intensity scaling
The band structure above is nice, but we can make the plot a little clearer by adjusting some of the parameters like 
the intensity scaling (via `--intensity`) and the colour map (via `--cmap`). Below we set `--intensity 2.5` to 
increase the colourmap intensity, and try `BuPu`, `viridis` and `bone_r` from left to right below: 

```bash
easyunfold unfold plot --intensity 2.5 --cmap "BuPu"
easyunfold unfold plot --intensity 2.5 --cmap "viridis"
easyunfold unfold plot --intensity 2.5 --cmap "bone_r"
```

BuPu             |  viridis        |  bone_r
:-------------------------:|:-------------------------:|:-------------------------:
![](../../examples/NaBiS2/NaBiS2_unfold-plot_BuPu.png)  |  ![](../../examples/NaBiS2/NaBiS2_unfold-plot_viridis.png) |  ![](../../examples/NaBiS2/NaBiS2_unfold-plot_bone_r.png)

## Unfolded Band Structure with Density of States (DOS)
We can plot the electronic density of states (DOS) alongside the unfolded band structure using the `--dos` option to 
provide the `vasprun.xml(.gz)` file from our supercell calculation:

```bash
easyunfold unfold plot --intensity 2 --dos vasprun.xml.gz --zero-line --dos-label DOS --gaussian 0.1
```

Here we've used some other plot options to customise the DOS plot; see the help message with 
`easyunfold unfold plot -h` for more info on this. 

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_dos.png
:alt: NaBiS2 unfolded band structure with DOS
:width: 700px

Unfolded band structure of NaBiS<sub>2</sub> alongside the electronic density of states (DOS)
```

:::{note}
For unfolded band structures calculated with (semi-)local DFT (LDA/GGA), you should not use the `vasprun.xml(.gz)` 
from the band structure calculation (which uses a non-uniform _k_-point mesh, thus giving an unrepresentative DOS 
output), but rather the preceding self-consistent calculation (used to obtain the `CHGCAR` for the LDA/GGA DOS 
calculation), or a separate DOS calculation.
:::

:::{tip}
To use the DOS plotting functionality of `easyunfold`, the `sumo` package must be installed. This is currently an 
optional dependency for `easyunfold` (to avoid strict requirements and `pymatgen` dependencies), but can be installed 
with `pip install sumo`.
:::

## Atom-Projected Unfolded Band Structure
We can also plot the unfolded band structure with atomic projections, which is useful for understanding the electronic 
structure of the material. In this case, we are curious as to which atoms are contributing to the band edges, and so 
the atomic projections will be useful. For this, we need the `PROCAR(.gz)` output from VASP with the atomic and orbital 
projection information, and so `LORBIT` should be set to `11`, `12`, `13` or `14` in the `INCAR` for the supercell 
calculation.

When plotting the unfolded band, the `plot-projections` subcommand is used with the `--procar <PROCAR>` and 
`--atoms` options:

```bash
easyunfold unfold plot-projections --atoms="Na,Bi,S" --intensity 3 --combined
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj.png
:alt: NaBiS2 atom-projected unfolded band structure
:width: 700px

Unfolded band structure of NaBiS<sub>2</sub> with atomic contributions.
```

:::{tip}
If the _k_-points have been split into multiple calculations (e.g. hybrid DFT band structures), the `--procar` option 
should be passed multiple times to specify the path to each split `PROCAR(.gz)` file (i.e. 
`--procar calc1/PROCAR --procar calc2/PROCAR ...`).
:::

From this plot, we can see that sulfur anions dominate the valence band, while bismuth cations dominate the conduction 
band, with minimal contributions from the sodium cations as expected.

:::{note}
The atomic projections are not stored in the `easyunfold.json` data file, so the `PROCAR(.gz)` file 
should be kept for replotting in the future.
:::

While the main conclusions of S dominating the valence band and Bi dominating the conduction band are clear from the 
plot above, the high intensity of the S projections in the valence band makes the Bi contributions in the conduction 
band more faint and less apparent. 

So, we can create separated plots for each atom type to make their individual contributions more clear, by omitting the 
`--combined` tag:

```bash
easyunfold unfold plot-projections --atoms="Na,Bi,S" --cmap="bone_r" --intensity 2
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_sep.png
:alt: NaBiS2 atom-projected separated unfolded band structure 
:width: 800px

Unfolded band structure of NaBiS<sub>2</sub> with atomic contributions plotted separately.
```

An alternative option here is also to just plot only the contributions of `Na` and `Bi` cations, with no S projections:
```bash
easyunfold unfold plot-projections --atoms="Na,Bi" --intensity 3 --combined
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_noS.png
:alt: NaBiS2 atom-projected unfolded band structure, no S 
:width: 700px

Unfolded band structure of NaBiS<sub>2</sub> with atomic contributions of only Na and Bi.
```

While this plot isn't the most aesthetic, it clearly shows that Bi (green) contributes to both the conduction band and 
(less so) valence band states, but Na (red) doesn't contribute significantly at or near the band edges 
(it's a spectator ion!). 

### Atom-projected Unfolded Band Structure with DOS
We can also combine the atom projections with the DOS plotting, using the `--dos` option as before:
```bash
easyunfold unfold plot-projections --atoms "Na,Bi,S" --intensity 3 --combined --dos vasprun.xml.gz --zero-line \
  --dos-label "DOS" --gaussian 0.1 --no-total --scale 2
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_dos.png
:alt: NaBiS2 atom-projected unfolded band structure with DOS
:width: 700px

Atom-projected unfolded band structure of NaBiS<sub>2</sub> alongside the electronic density of states (DOS)
```
The orbital contributions of elements in the DOS plot are automatically coloured to match that of the atomic 
projections in the unfolded band structure plot, and these colours can be changed with the `--colours` option (as shown 
in the [MgO example](https://smtg-Bham.github.io/easyunfold/examples/example_mgo.html)). 

## Unfolded Band Structure with Specific Atom Selection
In certain cases, we may want to project the contributions of specific atoms to the unfolded band structure, rather
than all atoms of a certain element. For example, in this example case of NaBiS$_2$, we see that a flat localised state 
is forming within the 'bulk' band gap. We can show that this state is contributed by only one or two atoms in the 
supercell, using the `--atoms-idx` option:
```bash
easyunfold unfold plot-projections -unfold plot-projections --atoms-idx="1-55,57-65,67-80|56,66" --intensity 2
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_sep_atomsidx.png
:alt: NaBiS2 atom-projected unfolded band structure, specific atoms 
:width: 700px

Unfolded band structure of NaBiS<sub>2</sub> with the atomic contributions of atoms 56 and 66 plotted separately.
```
Here we can see that this in-gap localised state is dominated by only two (sulfur) atoms in the supercell (atoms `56` 
and `66`), which correspond to Na-rich pockets in NaBiS<sub>2</sub>, as discussed in Huang & Kavanagh et al. [^1].

:::{note}
The `--atoms-idx` option is used to specify the atoms to be projected onto the unfolded band structure. This takes a 
string of the form `a-b|c-d|e-f` where `a`, `b`, `c`, `d`, `e` and `f` are integers corresponding to the atom indices 
in the VASP structure file (i.e. `POSCAR`/`CONTCAR`, corresponding to the `PROCAR(.gz)` being used to obtain the 
projections). Different groups are separated by `|`, and `-` can be used to define the range for each projected atom 
type. A comma-separated list can also be used instead of ranges with hyphens. Note that 1-based indexing is used for 
atoms, matching the convention in VASP, which is then converted to zero-based indexing internally in python. 
:::

## Unfolded Band Structure with Orbital Projections
If we want to see the contributions of specific orbitals to the unfolded band structure, we can use the `--orbitals`
option. This takes a string of the form `a,b|c|d,e,f` where `a`, `b`, `c`, `d`, `e` and `f` are the orbitals we want
to project onto the unfolded band structure. 

For example, if we want to see the contributions of the Bi $s$, $p$ and S $s$ orbitals to the unfolded band structure, 
we can use the following command:

```bash
easyunfold unfold plot-projections --atoms "Bi,Bi,S" --orbitals="s|p|s"  --intensity 3  --combined \
  --dos vasprun.xml.gz --zero-line --dos-label "DOS" --gaussian 0.1 --no-total --scale 5
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_sps_dos.png
:alt: NaBiS2 atom-projected unfolded band structure with DOS and Bi s/p separated
:width: 700px

Orbital-projected unfolded band structure of NaBiS<sub>2</sub> alongside the electronic density of states (DOS)
```
Here we have separated out the contributions of Bi $s$ orbitals (which have some weak anti-bonding contributions to the upper 
'bulk' VBM as expected, due to the occupied Bi lone-pair; see refs [^1] and [^2]) and $p$ orbitals (which contribute to 
the conduction band and lower valence band). We have also shown only the $s$ orbital contributions of sulfur, which we
can see have minimal contributions to the electronic structure near the band gap.


### $lm$-decomposed Orbital Projections
We can also use the `--orbitals` option to project onto specific $lm$-decomposed orbitals, and/or use the 
`--dos-orbitals` option to split the DOS plot into $lm$-decomposed orbital contributions. For example, if we want to
see the contributions of the Bi and S $p_x$, $p_y$ and $p_z$ orbitals to the unfolded band structure, we can use the 
following command:

```bash
easyunfold unfold plot-projections --atoms "Na,Bi,S" --orbitals="all|px,py,pz|px,py,pz" --intensity 3 --combined \
  --dos vasprun.xml.gz --zero-line --dos-label "DOS" --gaussian 0.1 --no-total --scale 6
```

```{figure} ../../examples/NaBiS2/NaBiS2_unfold-plot_proj_orbital_lm_dos.png
:alt: NaBiS2 $lm$-decomposed orbtial-projected unfolded band structure with DOS
:width: 700px

$lm$-decomposed orbital-projected unfolded band structure of NaBiS<sub>2</sub> alongside the DOS
```

Here, the $p_x$, $p_y$ and $p_z$ orbitals do not differ significantly in their contributions to the electronic 
structure, due to the cubic symmetry of the NaBiS<sub>2</sub> crystal structure. However, in other cases, such as
for the $d$ orbitals of transition metals in octahedral/tetrahedral environments, we would expect to see significant 
differences in the contributions of different $lm$-decomposed orbitals to the electronic structure.

:::{tip}
There are _many_ customisation options available for the plotting functions in `easyunfold`. See `easyunfold plot -h` or 
`easyunfold unfold plot-projections -h` for more details!
:::

[^1]: [Huang, YT., Kavanagh, S.R., Righetto, M. et al. Strong absorption and ultrafast localisation in NaBiS2 nanocrystals with slow charge-carrier recombination. Nat Commun 13, 4960 (2022)](https://www.nature.com/articles/s41467-022-32669-3) 
[^2]: [Wang, Y., Kavanagh, S.R., Burgués-Ceballos, I. et al. Cation disorder engineering yields AgBiS2 nanocrystals with enhanced optical absorption for efficient ultrathin solar cells. Nat. Photon. 16, 235–241 (2022).](https://www.nature.com/articles/s41566-021-00950-4)