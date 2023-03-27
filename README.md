# Easyunfold


![build](https://github.com/SMTG-UCL/easyunfold/actions/workflows/ci.yaml/badge.svg)
[![docs](https://github.com/SMTG-UCL/easyunfold/actions/workflows/docs.yaml/badge.svg)](https://smtg-ucl.github.io/easyunfold/)
[![codecov](https://codecov.io/gh/SMTG-UCL/easyunfold/branch/dev/graph/badge.svg?token=XLLWWU5UM2)](https://codecov.io/gh/SMTG-UCL/easyunfold)

![easyunfold](docs/img/logo.svg)

Documentation: https://smtg-ucl.github.io/easyunfold/

This package is intended for obtaining the effective band structure of a supercell for a certain path of the primitive cell.
It was originally based on [PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding) for reading wavefunction output of VASP, and contains some code of the latter.
An notable improvement is that breaking symmetry is taken accounted of by sampling additional kpoints and taking average accordingly, which was previously missing.
Our goal is to make the unfolding process easier to carry out and less likely to go wrong.

For the methodology, see: https://link.aps.org/doi/10.1103/PhysRevB.85.085201


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

## Installation

### Install from pip

The package can be installed from `pip`

```
pip install easyunfold
```

This will also install the dependencies, if they are missing.

After installation, run `easyunfold` should give the following output:

```
Usage: easyunfold [OPTIONS] COMMAND [ARGS]...

  Tool for performing band unfolding

Options:
  --help  Show this message and exit.

Commands:
  generate  Generate the kpoints for sampling the supercell
  unfold    Perform unfolding and plotting
```

### Install from source
A recently version of `pip` is needed to do this, due to the use of new style `pyproject.toml` configuration file.
To upgrade your `pip`, do:

```
pip install -U pip
```

Assuming the package is in the `easyunfold` folder, use the following command to install:

```
pip install ./easyunfold
```

## Studies using `easyunfold`
We'll add papers that use `easyunfold` to this list as they come out!
- Y. T. Huang & S. R. Kavanagh et al. [_Nature Communications_](https://www.nature.com/articles/s41467-022-32669-3) 2022
- A. T. J. Nicolson et al. [_ChemRxiv_](https://chemrxiv.org/engage/chemrxiv/article-details/63a5d1ffa53ea69e935559e2) 2023
- Y. Wang & S. R. Kavanagh et al. [_Nature Photonics_](https://www.nature.com/articles/s41566-021-00950-4) 2022 (early version)

<!-- Others? -->


## Contributors

Code Contributors:
  [Bonan Zhu](https://github.com/zhubonan)
  [Seán Kavanagh](https://github.com/kavanase)
  [Adair Nicolson](https://github.com/https://github.com/adair-nicolson)

And those who helped in the development:

  [Joe Willis](https://github.com/joebesity)
  [David O. Scanlon](http://davidscanlon.com/?page_id=5)
