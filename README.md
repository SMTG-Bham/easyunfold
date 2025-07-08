![build](https://github.com/SMTG-Bham/easyunfold/actions/workflows/ci.yaml/badge.svg)
[![docs](https://github.com/SMTG-Bham/easyunfold/actions/workflows/docs.yaml/badge.svg)](https://smtg-Bham.github.io/easyunfold/)
[![codecov](https://codecov.io/gh/SMTG-Bham/easyunfold/branch/main/graph/badge.svg?token=XLLWWU5UM2)](https://codecov.io/gh/SMTG-Bham/easyunfold)
[![PyPI](https://img.shields.io/pypi/v/easyunfold)](https://pypi.org/project/easyunfold)
[![Downloads](https://img.shields.io/pypi/dm/easyunfold)](https://smtg-Bham.github.io/easyunfold/)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.05974/status.svg)](https://doi.org/10.21105/joss.05974)

[![easyunfold](docs/img/logo.svg)](https://smtg-Bham.github.io/easyunfold/)

`easyunfold` is intended for obtaining the effective band structure of a supercell for a certain _k_-point
path of the primitive cell. It was originally based on
[PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding) for reading VASP wavefunction outputs,
with a notable improvement being that symmetry-breaking is properly accounted for by sampling necessary
additional _k_-points and averaging accordingly. Documentation site
[here](https://smtg-bham.github.io/easyunfold/)!

Our goal is to implement the band structure unfolding workflow in a robust and user-friendly software
package.
Typical applications of band structure unfolding are the electronic structure analysis of defects, disorder, alloys, surfaces (and more), as illustrated in the example outputs below and [docs examples](https://smtg-bham.github.io/easyunfold/examples.html).

For the methodology of supercell band unfolding, see
[here](https://link.aps.org/doi/10.1103/PhysRevB.85.085201).

### Example Outputs
[Cs₂(Sn/Ti)Br₆ Vacancy-Ordered Perovskite Alloys](https://doi.org/10.1021/acs.jpcc.3c05204) |     Oxygen Vacancy (*V*ₒ⁰) in MgO
:-------------------------:|:------------------------------------:
<img src="docs/img/CSTB_easyunfold.gif" height="400"/> | <img src="examples/MgO/unfold_project_MgO_v_O_0_tall.png" height="400"/>

See the [`easyunfold` YouTube tutorial](https://youtu.be/9zeABbd1r1U?si=Oix3Bamiw8DZaMO4) for a quick overview of the theory of band structure unfolding, and a walkthrough of the calculation & analysis workflow with `easyunfold`.

## Usage

To generate an unfolded band structure, one typically needs to perform the following steps:

1. Create a primitive unit cell, and generate a band structure _k_-point path corresponding to this
   primitive cell.
2. Create a supercell (e.g. disordered, defective, surface slab etc.), and obtain its optimised structure.
3. Generate a series of _k_-points for the supercell to be calculated.
4. Perform a band structure calculation with the supercell, and save its wavefunction output to file.
5. Post-process the supercell wavefunction to obtain the unfolded band structure in the _k_-point path
   of the primitive unit cell.

These generation and analysis steps are automated in `easyunfold`, with only the primitive unit cell and
supercell structures required as inputs from the user.

Typically, the supercell comprises some form of symmetry-breaking relative to the primitive cell, such
as defects, disorder (e.g. special quasi-random structures (SQS) for site disorder – other forms of
disorder such as magnetic, dynamic/vibrational, polar, elastic etc. also possible), or a surface/interface
slab.
In all cases, the supercell symmetry is lowered compared to the pristine primitive cell.
Hence, for a given _k_-point path in the primitive cell Brillouin Zone, additional _k_-points are
required to be sampled for the supercell, and the extracted spectral weights need to be appropriately
averaged to obtain the correct effective band structure (EBS). See the docs
[Theory](https://smtg-bham.github.io/easyunfold/theory.html) page for more details.
<!-- when JOSS submitted, add link to paper (discussion of theory) here! -->
<!--- When JOSS submitted, add 'License and Citation' section here, and `CITATION.cff` file --->

Please see the [documentation](https://smtg-bham.github.io/easyunfold/) for guides and examples.

## Installation

### Install from `pip`

`easyunfold` can be installed from `pip`:

```
pip install easyunfold
```

This will also install the package dependencies, if any are missing.

After installation, running `easyunfold` on the command-line should give the following output:

```
Usage: easyunfold [OPTIONS] COMMAND [ARGS]...

  Tool for performing band unfolding

Options:
  --help  Show this message and exit.

Commands:
  generate  Generate the kpoints for sampling the supercell
  unfold    Perform unfolding and plotting
```

### Developer Installation (from source)

A recent version of `pip` is needed to do this, due to the new style of the `pyproject.toml` configuration
file.
To upgrade your `pip`, do:

```
pip install -U pip
```

Assuming the package is in the `easyunfold` folder, use the following command to install:

```
pip install "./easyunfold[test,doc,pre-commit]"
```

which also installs additional dependencies for building documentation (`doc`), running tests (`test`) and
dependencies for using pre-commit hooks (`pre-commit`).

## Citation

If you use `easyunfold` in your work, please cite:
- B. Zhu, S. R. Kavanagh & D. O. Scanlon, (2024). easyunfold: A Python package for unfolding electronic band structures. Journal of Open Source Software, 9(93), 5974, https://doi.org/10.21105/joss.05974

## Studies using `easyunfold`

We'll add papers that use `easyunfold` to this list as they come out!

- S. Husremović et al. **_Local interface effects modulate global charge order
and optical properties of 1T-TaS<sub>2</sub>/1H-WSe<sub>2</sub> heterostructures_** [_arXiv_](https://arxiv.org/abs/2508.01512) 2025
- L. Richarz et al. **_Ferroelectric domain walls for environmental sensors_** [_ACS Applied Materials & Interfaces_](http://dx.doi.org/10.1021/acsami.5c04875) 2025
- J. M. Domínguez-Vázquez et al. **_Thermoelectric performance boost by chemical order in epitaxial L2<sub>1</sub> (100) and (110) oriented undoped Fe<sub>2</sub>VAl thin films: an experimental and theoretical study_** [_Journal of Materials Chemistry A_](https://doi.org/10.1039/D5TA02619A) 2025
- J. Tu et al. **_Giant switchable ferroelectric photovoltage in double-perovskite epitaxial films through chemical negative strain_** [_Science Advances_](https://doi.org/10.1126/sciadv.ads4925) 2025
- L. F. Leon-Pinzon et al. **_Observation of Pseudogap in Cr<sub>1−x</sub>Y<sub>x</sub>N magnetic alloy and its impact on the Seebeck coefficient by ab-initio calculations_** [_arXiv_](http://dx.doi.org/10.48550/arXiv.2506.00687) 2025
- C. Zhang & J. Recatala-Gomez & Z. Aabdin & Y. Jiang et al. **_Direct Joule-Heated Non-Equilibrium Synthesis Enables High Performing Thermoelectrics_** [_arXiv_](https://arxiv.org/abs/2506.04447) 2025
- W. Feng et al. **_Unraveling the role of oxygen vacancies in the electronic and optical properties of κ-Ga2O3_** [_ResearchSquare_](https://doi.org/10.21203/rs.3.rs-6624989/v1) 2025
- J. J. Plata et al. **_High-Entropy Skutterudites as Thermoelectrics: Synthesizability and Band Convergence via the Cocktail Effect_** [_ChemRxiv_](https://chemrxiv.org/engage/chemrxiv/article-details/67ff672750018ac7c59b7b11) 2025
- P. B. Colominas et al. **_Giant Thermally Induced Band-Gap Renormalization in Anharmonic Silver Chalcohalide Antiperovskites_** [_Journal of Materials Chemistry C_](https://doi.org/10.1039/D5TC00863H) 2025
- L. Zhang et al. **_Mg doping point defects in Al<sub>0.5</sub>Ga<sub>0.5</sub>N by density functional theory_** [_Vacuum_](https://doi.org/10.1016/j.vacuum.2025.114610) 2025
- L. Zhang et al. **_Impurity point defects in Mg doping Al<sub>0.5</sub>Ga<sub>0.5</sub>N: A first principles study_** [_Computational Materials Science_](https://doi.org/10.1016/j.commatsci.2025.113925) 2025
- L. Zhang et al. **_Study of native point defects in Al<sub>0.5</sub>Ga<sub>0.5</sub>N by first principles calculations_** [_Computational Materials Science_](https://doi.org/10.1016/j.commatsci.2024.113312) 2024
- H. Maleki-Ghaleh et al. **_Visible Light-Sensitive Sustainable Quantum Dot Crystals of Co/Mg Doped Natural Hydroxyapatite Possessing Antimicrobial Activity and Biocompatibility_** [_Small_](https://doi.org/10.1002/smll.202405708) 2024
- K. Eggestad, B. A. D. Williamson, D. Meier and S. M. Selbach **_Mobile Intrinsic Point Defects for Conductive Neutral Domain Walls in LiNbO<sub>3</sub>_** [_Journal of Materials Chemistry C_](https://doi.org/10.1039/D4TC02856B) 2024
- Dargahi et al. **_Synthesis and characterization of zinc-doped hematite nanoparticles for photocatalytic applications and their electronic structure studies by density functional theory_** [_Optical Materials_](https://doi.org/10.1016/j.optmat.2024.116234) 2024
- S. M. Liga & S. R. Kavanagh, A. Walsh, D. O. Scanlon and G. Konstantatos **_Mixed-Cation Vacancy-Ordered Perovskites (Cs<sub>2</sub>Ti<sub>1–x</sub>Sn<sub>x</sub>X<sub>6</sub>; X = I or Br): Low-Temperature Miscibility, Additivity, and Tunable Stability_** [_Journal of Physical Chemistry C_](https://doi.org/10.1021/acs.jpcc.3c05204) 2023
- Y. T. Huang & S. R. Kavanagh et al. **_Strong absorption and ultrafast localisation in NaBiS<sub>2</sub> nanocrystals with slow charge-carrier recombination_** [_Nature Communications_](https://www.nature.com/articles/s41467-022-32669-3) 2022
- A. T. J. Nicolson et al. **_Interplay of Static and Dynamic Disorder in the Mixed-Metal Chalcohalide Sn<sub>2</sub>SbS<sub>2</sub>I<sub>3</sub>_** [_Journal of the Americal Chemical Society_](https://doi.org/10.1021/jacs.2c13336) 2023
- Y. Wang & S. R. Kavanagh et al. **_Cation disorder engineering yields AgBiS<sub>2</sub> nanocrystals with enhanced optical absorption for efficient ultrathin solar cells_** [_Nature Photonics_](https://www.nature.com/articles/s41566-021-00950-4) 2022 (early version)
<!-- Others? -->

## DFT code support

At the moment, easyunfold supports [VASP](https://www.vasp.at) and [CASTEP](http://www.castep.org), but most of the routines are abstracted from the code specific details.
In principle, support for other plane wave DFT code can be added by:

- Implementing a subclass of `WaveFunction` that handles reading the wave function output.
- Implementing functions for reading/writing _k_-points.
- Adding branches for dispatching based on the `dft_code` attribute of the `UnfoldKSet` object in
  various places within the code.

The Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) is used by `easyunfold` for
reading in structures, so structure file IO is natively supported for essentially all public DFT codes.

### Code Compatibility Notes
- Atom-projected band structures are currently only supported for `VASP` calculation outputs.
- Gamma-only and non-collinear spin calculations are not supported for `CASTEP`.

## Contributors
- [Bonan Zhu](https://github.com/zhubonan)
- [Seán Kavanagh](https://github.com/kavanase)
- [Adair Nicolson](https://github.com/adair-nicolson)

And those who helped in the development:
- [Joe Willis](https://github.com/joebesity)
- [David O. Scanlon](http://davidscanlon.com/?page_id=5)
