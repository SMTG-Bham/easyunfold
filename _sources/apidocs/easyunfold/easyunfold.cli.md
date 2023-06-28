# {py:mod}`easyunfold.cli`

```{py:module} easyunfold.cli
```

```{autodoc2-docstring} easyunfold.cli
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`easyunfold <easyunfold.cli.easyunfold>`
  - ```{autodoc2-docstring} easyunfold.cli.easyunfold
    :summary:
    ```
* - {py:obj}`generate <easyunfold.cli.generate>`
  - ```{autodoc2-docstring} easyunfold.cli.generate
    :summary:
    ```
* - {py:obj}`unfold <easyunfold.cli.unfold>`
  - ```{autodoc2-docstring} easyunfold.cli.unfold
    :summary:
    ```
* - {py:obj}`unfold_status <easyunfold.cli.unfold_status>`
  - ```{autodoc2-docstring} easyunfold.cli.unfold_status
    :summary:
    ```
* - {py:obj}`unfold_calculate <easyunfold.cli.unfold_calculate>`
  - ```{autodoc2-docstring} easyunfold.cli.unfold_calculate
    :summary:
    ```
* - {py:obj}`add_plot_options <easyunfold.cli.add_plot_options>`
  - ```{autodoc2-docstring} easyunfold.cli.add_plot_options
    :summary:
    ```
* - {py:obj}`unfold_effective_mass <easyunfold.cli.unfold_effective_mass>`
  - ```{autodoc2-docstring} easyunfold.cli.unfold_effective_mass
    :summary:
    ```
* - {py:obj}`unfold_plot <easyunfold.cli.unfold_plot>`
  - ```{autodoc2-docstring} easyunfold.cli.unfold_plot
    :summary:
    ```
* - {py:obj}`unfold_plot_projections <easyunfold.cli.unfold_plot_projections>`
  - ```{autodoc2-docstring} easyunfold.cli.unfold_plot_projections
    :summary:
    ```
* - {py:obj}`_unfold_plot <easyunfold.cli._unfold_plot>`
  - ```{autodoc2-docstring} easyunfold.cli._unfold_plot
    :summary:
    ```
* - {py:obj}`print_symmetry_data <easyunfold.cli.print_symmetry_data>`
  - ```{autodoc2-docstring} easyunfold.cli.print_symmetry_data
    :summary:
    ```
* - {py:obj}`matrix_from_string <easyunfold.cli.matrix_from_string>`
  - ```{autodoc2-docstring} easyunfold.cli.matrix_from_string
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`SUPPORTED_DFT_CODES <easyunfold.cli.SUPPORTED_DFT_CODES>`
  - ```{autodoc2-docstring} easyunfold.cli.SUPPORTED_DFT_CODES
    :summary:
    ```
* - {py:obj}`DEFAULT_CMAPS <easyunfold.cli.DEFAULT_CMAPS>`
  - ```{autodoc2-docstring} easyunfold.cli.DEFAULT_CMAPS
    :summary:
    ```
* - {py:obj}`CONTEXT_SETTINGS <easyunfold.cli.CONTEXT_SETTINGS>`
  - ```{autodoc2-docstring} easyunfold.cli.CONTEXT_SETTINGS
    :summary:
    ```
````

### API

````{py:data} SUPPORTED_DFT_CODES
:canonical: easyunfold.cli.SUPPORTED_DFT_CODES
:value: >
   ('vasp', 'castep')

```{autodoc2-docstring} easyunfold.cli.SUPPORTED_DFT_CODES
```

````

````{py:data} DEFAULT_CMAPS
:canonical: easyunfold.cli.DEFAULT_CMAPS
:value: >
   ['Purples', 'Greens', 'Oranges', 'Reds', 'Blue', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',...

```{autodoc2-docstring} easyunfold.cli.DEFAULT_CMAPS
```

````

````{py:data} CONTEXT_SETTINGS
:canonical: easyunfold.cli.CONTEXT_SETTINGS
:value: >
   None

```{autodoc2-docstring} easyunfold.cli.CONTEXT_SETTINGS
```

````

````{py:function} easyunfold()
:canonical: easyunfold.cli.easyunfold

```{autodoc2-docstring} easyunfold.cli.easyunfold
```
````

````{py:function} generate(pc_file, code, sc_file, matrix, kpoints, time_reversal, out_file, no_expand, symprec, nk_per_split, scf_kpoints, yes)
:canonical: easyunfold.cli.generate

```{autodoc2-docstring} easyunfold.cli.generate
```
````

````{py:function} unfold(ctx, data_file, mpl_style_file)
:canonical: easyunfold.cli.unfold

```{autodoc2-docstring} easyunfold.cli.unfold
```
````

````{py:function} unfold_status(ctx)
:canonical: easyunfold.cli.unfold_status

```{autodoc2-docstring} easyunfold.cli.unfold_status
```
````

````{py:function} unfold_calculate(ctx, wavefunc, save_as, gamma, ncl)
:canonical: easyunfold.cli.unfold_calculate

```{autodoc2-docstring} easyunfold.cli.unfold_calculate
```
````

````{py:function} add_plot_options(func)
:canonical: easyunfold.cli.add_plot_options

```{autodoc2-docstring} easyunfold.cli.add_plot_options
```
````

````{py:function} unfold_effective_mass(ctx, intensity_threshold, spin, band_filter, npoints, extrema_detect_tol, degeneracy_detect_tol, nocc, plot, plot_fit, fit_label, out_file)
:canonical: easyunfold.cli.unfold_effective_mass

```{autodoc2-docstring} easyunfold.cli.unfold_effective_mass
```
````

````{py:function} unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar, atoms_idx, orbitals, title, width, height, dpi)
:canonical: easyunfold.cli.unfold_plot

```{autodoc2-docstring} easyunfold.cli.unfold_plot
```
````

````{py:function} unfold_plot_projections(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar, atoms_idx, orbitals, title, combined, intensity, colors, width, height, dpi)
:canonical: easyunfold.cli.unfold_plot_projections

```{autodoc2-docstring} easyunfold.cli.unfold_plot_projections
```
````

````{py:function} _unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar, atoms_idx, orbitals, title, width, height, dpi, ax=None)
:canonical: easyunfold.cli._unfold_plot

```{autodoc2-docstring} easyunfold.cli._unfold_plot
```
````

````{py:function} print_symmetry_data(kset)
:canonical: easyunfold.cli.print_symmetry_data

```{autodoc2-docstring} easyunfold.cli.print_symmetry_data
```
````

````{py:function} matrix_from_string(string)
:canonical: easyunfold.cli.matrix_from_string

```{autodoc2-docstring} easyunfold.cli.matrix_from_string
```
````
