# Easyunfold documentation

![build](https://github.com/SMTG-UCL/easyunfold/actions/workflows/ci.yaml/badge.svg)
[![docs](https://github.com/SMTG-UCL/easyunfold/actions/workflows/docs.yaml/badge.svg)](https://smtg-ucl.github.io/easyunfold/)
[![codecov](https://codecov.io/gh/SMTG-UCL/easyunfold/branch/main/graph/badge.svg?token=XLLWWU5UM2)](https://codecov.io/gh/SMTG-UCL/easyunfold)

This package is intended for obtaining the effective band structure of a supercell for a certain path of the primitive cell.
It was originally based on [PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding) for reading wavefunction output of VASP, and contains some code of the latter.
An notable improvement is that breaking symmetry is taken accounted of by sampling additional kpoints and taking average accordingly, which was previously missing.
Our goal is to make the unfolding process easier to carry out and less likely to go wrong.

For the methodology, see [here](https://link.aps.org/doi/10.1103/PhysRevB.85.085201).


## Usage

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

At the moment, only VASP calculations are supported, although in principle other PW code can be supported easily if the wavefunction can be read in.

Please see the documentation for guide and examples.

## Studies using `easyunfold`

We'll add papers that use `easyunfold` to this list as they come out!

- Y. T. Huang & S. R. Kavanagh et al. [_Nature Communications_](https://www.nature.com/articles/s41467-022-32669-3) 2022
- A. T. J. Nicolson et al. [_ChemRxiv_](https://chemrxiv.org/engage/chemrxiv/article-details/63a5d1ffa53ea69e935559e2) 2023
- Y. Wang & S. R. Kavanagh et al. [_Nature Photonics_](https://www.nature.com/articles/s41566-021-00950-4) 2022 (early version)

## DFT code support

At the moment, easyunfold supports VASP and CASTEP, but most of the routines are abstracted from the code specific details.
In principle, support for other plane wave DFT code can be added by:

- Implementing a subclass of `WaveFunction` that handles reading the wave function output.
- Implementing functions for reading/writing kpoints.
- Adding branches for dispatching based on the `dft_code` attribute of the `UnfoldKSet` object in various places.

The Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) is used by this package for reading in structures, so there is usually need to implement specialised parser for structure file.


## Known limitations

- Atomic projection is available for VASP calculation only for now.
- CASTEP gamma-only and non-collinear spin calculations are not supported. 

## Contributors

Code Contributors:
  [Bonan Zhu](https://github.com/zhubonan)  
  [Seán Kavanagh](https://github.com/kavanase)  
  [Adair Nicolson](https://github.com/adair-nicolson)  

And those who helped in the development:

  [Joe Willis](https://github.com/joebesity)  
  [David O. Scanlon](http://davidscanlon.com/?page_id=5)  

## Contributing


## Bugs reports and feature requests


Bug reports and feature requests are well come.
If you found any bug or missing features please report it on the [Issue Tracker](https://github.com/SMTG-UCL/easyunfold/issues).
Poor or incorrect documentation is considered a bug, but please be as specific as
possible when asking for improvements.

## Code contributions

We welcome your help in improving and extending the package with your
own contributions. This is managed through Github pull requests;
for external contributions, we prefer the
[fork and pull](https://guides.github.com/activities/forking/)
workflow while core developers use branches in the main repository:

1. First open an Issue to discuss the proposed contribution. This
  discussion might include how the changes fit easyunfold's scope and a
  general technical approach.
2. Make your own project fork and implement the changes
  there. Please keep your code style compliant and use the pre-commit hooks.
3. Open a pull request to merge the changes into the main
  project. A more detailed discussion can take place there before
  the changes are accepted.


# Table of contents

```{toctree}
:maxdepth: 2

installation.md
guide.md
examples.md
theory.md
cli.rst
apidocs/index
```