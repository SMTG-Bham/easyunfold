# {py:mod}`easyunfold.plotting`

```{py:module} easyunfold.plotting
```

```{autodoc2-docstring} easyunfold.plotting
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`UnfoldPlotter <easyunfold.plotting.UnfoldPlotter>`
  - ```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`interpolate_colors <easyunfold.plotting.interpolate_colors>`
  - ```{autodoc2-docstring} easyunfold.plotting.interpolate_colors
    :summary:
    ```
````

### API

`````{py:class} UnfoldPlotter(unfold: easyunfold.unfold.UnfoldKSet)
:canonical: easyunfold.plotting.UnfoldPlotter

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter.__init__
```

````{py:method} plot_spectral_function(engs: numpy.ndarray, sf: numpy.ndarray, eref=None, figsize=(4, 3), ylim=(-3, 3), dpi=150, vscale=1.0, contour_plot=False, alpha=1.0, save=False, ax=None, vmin=None, vmax=None, cmap='PuRd', show=False, title=None)
:canonical: easyunfold.plotting.UnfoldPlotter.plot_spectral_function

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter.plot_spectral_function
```

````

````{py:method} _plot_spectral_function_rgba(engs: numpy.ndarray, sf: numpy.ndarray, eref: typing.Union[None, float] = None, figsize=(4, 3), ylim=(-3, 3), dpi: float = 150, intensity: float = 1.0, save: bool = False, ax: typing.Union[None, matplotlib.pyplot.Axes] = None, show: bool = False, title: typing.Union[None, str] = None, vmin: typing.Union[None, float] = None, vmax: typing.Union[None, float] = None)
:canonical: easyunfold.plotting.UnfoldPlotter._plot_spectral_function_rgba

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter._plot_spectral_function_rgba
```

````

````{py:method} _add_kpoint_labels(ax: matplotlib.pyplot.Axes, x_is_kidx=False)
:canonical: easyunfold.plotting.UnfoldPlotter._add_kpoint_labels

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter._add_kpoint_labels
```

````

````{py:method} plot_effective_mass(eff: easyunfold.effective_mass.EffectiveMass, engs: numpy.ndarray, sf: numpy.ndarray, eref: typing.Union[None, float] = None, save: typing.Union[None, str] = None, show: bool = False, effective_mass_data: dict = None, **kwargs)
:canonical: easyunfold.plotting.UnfoldPlotter.plot_effective_mass

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter.plot_effective_mass
```

````

````{py:method} plot_spectral_weights(figsize=(4, 3), ylim=(-3, 3), dpi: float = 150, factor: float = 3.0, eref: typing.Union[None, float] = None, color: str = 'C1', alpha: float = 0.5, save: typing.Union[None, str] = None, ax: typing.Union[None, matplotlib.pyplot.Axes] = None, show: bool = False, title: typing.Union[None, str] = None)
:canonical: easyunfold.plotting.UnfoldPlotter.plot_spectral_weights

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter.plot_spectral_weights
```

````

````{py:method} plot_projected(procar: typing.Union[str, list], eref=None, gamma=False, npoints=2000, sigma=0.2, ncl=False, symm_average=True, figsize=(4, 3), ylim=(-3, 3), dpi=150, vscale=1.0, contour_plot=False, alpha=1.0, save=False, ax=None, vmin=None, vmax=None, cmap='PuRd', show=False, title=None, atoms_idx=None, orbitals=None, intensity=1.0, use_subplot=False, colors=['r', 'g', 'b', 'purple'], colorspace='lab')
:canonical: easyunfold.plotting.UnfoldPlotter.plot_projected

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter.plot_projected
```

````

````{py:method} plot_effective_mass_fit(efm: easyunfold.effective_mass.EffectiveMass, npoints: int = 3, carrier: str = 'electrons', idx: int = 0, ax: typing.Union[matplotlib.pyplot.Axes, None] = None, save: typing.Union[None, str] = None, dpi: float = 120)
:canonical: easyunfold.plotting.UnfoldPlotter.plot_effective_mass_fit
:staticmethod:

```{autodoc2-docstring} easyunfold.plotting.UnfoldPlotter.plot_effective_mass_fit
```

````

`````

````{py:function} interpolate_colors(colors: list, weights: list, colorspace='lab', normalize=True)
:canonical: easyunfold.plotting.interpolate_colors

```{autodoc2-docstring} easyunfold.plotting.interpolate_colors
```
````
