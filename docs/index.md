# `easyunfold` Documentation

![build](https://github.com/SMTG-Bham/easyunfold/actions/workflows/ci.yaml/badge.svg)
[![docs](https://github.com/SMTG-Bham/easyunfold/actions/workflows/docs.yaml/badge.svg)](https://smtg-Bham.github.io/easyunfold/)
[![codecov](https://codecov.io/gh/SMTG-Bham/easyunfold/branch/main/graph/badge.svg?token=XLLWWU5UM2)](https://codecov.io/gh/SMTG-Bham/easyunfold)
[![PyPI](https://img.shields.io/pypi/v/easyunfold)](https://pypi.org/project/easyunfold)
[![Downloads](https://img.shields.io/pypi/dm/easyunfold)](https://smtg-Bham.github.io/easyunfold/)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.05974/status.svg)](https://doi.org/10.21105/joss.05974)

This package is intended for obtaining the effective band structure of a supercell for a certain _k_-point
path of the primitive cell. It was originally based on 
[PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding) for reading VASP wavefunction outputs, 
with a notable improvement being that symmetry-breaking is properly accounted for by sampling necessary 
additional _k_-points and averaging accordingly.
Typical applications of band structure unfolding are the electronic structure analysis of defects, disorder, alloys, surfaces (and more), as illustrated in the example outputs below and [docs examples](https://smtg-Bham.github.io/easyunfold/examples.html).

Our goal is to implement the band structure unfolding workflow in a robust and user-friendly software 
package.

For the methodology of supercell band unfolding, see 
[here](https://link.aps.org/doi/10.1103/PhysRevB.85.085201).

## Example Outputs
| [Cs₂(Sn/Ti)Br₆ Vacancy-Ordered Perovskite Alloys](https://doi.org/10.1021/acs.jpcc.3c05204) |                        Oxygen Vacancy (*V*ₒ⁰) in MgO                         |
|:-------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------:|
|                      <img src="img/CSTB_easyunfold.gif" height="400"/>                      | <img src="../examples/MgO/unfold_project_MgO_v_O_0_tall.png" height="400"/> |

:::{tip}
See the [`easyunfold` YouTube tutorial](https://youtu.be/9zeABbd1r1U?si=Oix3Bamiw8DZaMO4) for a quick overview of the theory of band structure unfolding, and a walkthrough of the calculation & analysis workflow with `easyunfold`. 
:::

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
[Theory](https://smtg-Bham.github.io/easyunfold/theory.html) page and/or [JOSS paper](https://doi.org/10.21105/joss.05974) for more details.

## Citation

If you use `easyunfold` in your work, please cite:
- B. Zhu, S. R. Kavanagh & D. O. Scanlon, (2024). easyunfold: A Python package for unfolding electronic band structures. Journal of Open Source Software, 9(93), 5974, https://doi.org/10.21105/joss.05974

## Studies using `easyunfold`

We'll add papers that use `easyunfold` to this list as they come out!

- P. Russell et al. **_Computational prediction of Y-doped Cd<sub>2</sub>Sb<sub>2</sub>O<sub>7</sub> as a competitive Sb-based n-type Transparent Conducting Oxide_** [_ChemRxiv_](https://doi.org/10.26434/chemrxiv-2025-c1r3l) 2025
- R. E. Philip et al. **_Disorder mediated fully compensated ferrimagnetic spin-gapless semiconducting behaviour in Cr<sub>3</sub>Al Heusler alloy_** [_arXiv_](https://arxiv.org/abs/2512.10885) 2025
- S. Husremović et al. **_Local interface effects modulate global charge order and optical properties of 1T-TaS<sub>2</sub>/1H-WSe<sub>2</sub> heterostructures_** [_arXiv_](https://arxiv.org/abs/2508.01512) 2025
- L. Richarz et al. **_Ferroelectric domain walls for environmental sensors_** [_ACS Applied Materials & Interfaces_](http://dx.doi.org/10.1021/acsami.5c04875) 2025
- J. M. Domínguez-Vázquez et al. **_Thermoelectric performance boost by chemical order in epitaxial L2<sub>1</sub> (100) and (110) oriented undoped Fe<sub>2</sub>VAl thin films: an experimental and theoretical study_** [_Journal of Materials Chemistry A_](https://doi.org/10.1039/D5TA02619A) 2025
- J. Tu et al. **_Giant switchable ferroelectric photovoltage in double-perovskite epitaxial films through chemical negative strain_** [_Science Advances_](https://doi.org/10.1126/sciadv.ads4925) 2025
- L. F. Leon-Pinzon et al. **_Observation of Pseudogap in Cr<sub>1−x</sub>Y<sub>x</sub>N magnetic alloy and its impact on the Seebeck coefficient by ab-initio calculations_** [_arXiv_](http://dx.doi.org/10.48550/arXiv.2506.00687) 2025
- C. Zhang & J. Recatala-Gomez & Z. Aabdin & Y. Jiang et al. **_Direct Joule-Heated Non-Equilibrium Synthesis Enables High Performing Thermoelectrics_** [_arXiv_](https://arxiv.org/abs/2506.04447) 2025
- W. Feng et al. **_Unraveling the role of oxygen vacancies in the electronic and optical properties of κ-Ga2O3_** [_ResearchSquare_](https://doi.org/10.21203/rs.3.rs-6624989/v1) 2025
- J. J. Plata et al. **_High-Entropy Skutterudites as Thermoelectrics: Synthesizability and Band Convergence via the Cocktail Effect_** [_ChemRxiv_](https://chemrxiv.org/engage/chemrxiv/article-details/67ff672750018ac7c59b7b11) 2025
- Y. Gong Wang et al. **_Influence of Vanadium and Chromium Doping on the Thermoelectric Performance of AgSbTe<sub>2</sub>_** [_Physica Scripta_](https://doi.org/10.1088/1402-4896/ae26ec) 2025
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

## DFT code support

At the moment, easyunfold supports [VASP](https://www.vasp.at) and [CASTEP](http://www.castep.org), but most of the routines are abstracted from the code specific details.
In principle, support for other plane wave DFT code can be added by:

- Implementing a subclass of `WaveFunction` that handles reading the wave function output.
- Implementing functions for reading/writing _k_-points.
- Adding branches for dispatching based on the `dft_code` attribute of the `UnfoldKSet` object in 
  various places within the code.

The Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) is used by `easyunfold` for 
reading in structures, so structure file IO is natively supported for essentially all public DFT codes.

In fact, ASE can already run band structure calculations using many plane-wave DFT codes.
However, reading the plane wave coefficients from calculation outputs is not widely supported yet, which are needed here for band unfolding.
Nevertheless, using ASE's existing IO framework to widen the code support can be a fruitful direction for further development. 

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


## Bugs reports and feature requests
Bug reports and feature requests are well come.
If you found any bug or missing features please report it on the 
[Issue Tracker](https://github.com/SMTG-Bham/easyunfold/issues).

## Seeking support 

If you need support about using the software, please open a ticket with the *help wanted* label on the [Issue Tracker](https://github.com/SMTG-Bham/easyunfold/issues).

## Contributing

### Code contributions
We welcome your help in improving and extending the package with your
own contributions. This is managed through GitHub pull requests;
for external contributions, we prefer the
[fork and pull](https://guides.github.com/activities/forking/)
workflow while core developers use branches in the main repository:

1. First open an [Issue](https://github.com/SMTG-Bham/easyunfold/issues) to discuss the proposed 
   contribution. This discussion might include how the changes fit `easyunfold`'s scope and a
  general technical approach.
2. Make your own project fork and implement the changes
  there. Please keep your code style compliant and use the `pre-commit` hooks.
3. Open a pull request to merge the changes into the main
  project. A more detailed discussion can take place there before
  the changes are accepted.

### Development 

To develope the package, please clone it on GitHub and install it with the `doc` and `test` extras so all dependencies needed for development are installed:

```bash
git clone https://github.com/SMTG-Bham/easyunfold
pip install -e "./easyunfold[doc,test]"
```

To run the tests, simply run the `pytest` command  while in the top-level code repository.

To build the documentation, run the `make html` command while in the `docs` directory.

If you want to add new features, please also include a test case to check that they work as expected.
For more information please consult [pytest's documentation](https://docs.pytest.org).

Please note that we use pre-commit hooks for code formatting and linting. 
Please install them using `pre-commit install` so these hooks can run before committing code updates.


```{toctree}
  :hidden:
  :maxdepth: 2

installation.md
tutorial.md
examples.md
theory.md
references.md
publications.md
cli.rst
apidocs/index
```