# Installation
## Install from `pip`

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
pip install ./easyunfold
```
