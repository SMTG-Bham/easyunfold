# {py:mod}`easyunfold.effective_mass`

```{py:module} easyunfold.effective_mass
```

```{autodoc2-docstring} easyunfold.effective_mass
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`EffectiveMass <easyunfold.effective_mass.EffectiveMass>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`fit_effective_mass <easyunfold.effective_mass.fit_effective_mass>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.fit_effective_mass
    :summary:
    ```
* - {py:obj}`fitted_band <easyunfold.effective_mass.fitted_band>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.fitted_band
    :summary:
    ```
* - {py:obj}`points_with_tol <easyunfold.effective_mass.points_with_tol>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.points_with_tol
    :summary:
    ```
* - {py:obj}`locate_kpoint_segment <easyunfold.effective_mass.locate_kpoint_segment>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.locate_kpoint_segment
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`eV_to_hartree <easyunfold.effective_mass.eV_to_hartree>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.eV_to_hartree
    :summary:
    ```
* - {py:obj}`bohr_to_m <easyunfold.effective_mass.bohr_to_m>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.bohr_to_m
    :summary:
    ```
* - {py:obj}`angstrom_to_bohr <easyunfold.effective_mass.angstrom_to_bohr>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.angstrom_to_bohr
    :summary:
    ```
* - {py:obj}`TMP_DATA <easyunfold.effective_mass.TMP_DATA>`
  - ```{autodoc2-docstring} easyunfold.effective_mass.TMP_DATA
    :summary:
    ```
````

### API

````{py:data} eV_to_hartree
:canonical: easyunfold.effective_mass.eV_to_hartree
:value: >
   None

```{autodoc2-docstring} easyunfold.effective_mass.eV_to_hartree
```

````

````{py:data} bohr_to_m
:canonical: easyunfold.effective_mass.bohr_to_m
:value: >
   None

```{autodoc2-docstring} easyunfold.effective_mass.bohr_to_m
```

````

````{py:data} angstrom_to_bohr
:canonical: easyunfold.effective_mass.angstrom_to_bohr
:value: >
   None

```{autodoc2-docstring} easyunfold.effective_mass.angstrom_to_bohr
```

````

````{py:data} TMP_DATA
:canonical: easyunfold.effective_mass.TMP_DATA
:value: >
   None

```{autodoc2-docstring} easyunfold.effective_mass.TMP_DATA
```

````

````{py:function} fit_effective_mass(distances, energies, parabolic=True)
:canonical: easyunfold.effective_mass.fit_effective_mass

```{autodoc2-docstring} easyunfold.effective_mass.fit_effective_mass
```
````

````{py:function} fitted_band(x: numpy.ndarray, eff_mass: float) -> numpy.ndarray
:canonical: easyunfold.effective_mass.fitted_band

```{autodoc2-docstring} easyunfold.effective_mass.fitted_band
```
````

````{py:function} points_with_tol(array, value, tol=0.0001, sign=1)
:canonical: easyunfold.effective_mass.points_with_tol

```{autodoc2-docstring} easyunfold.effective_mass.points_with_tol
```
````

`````{py:class} EffectiveMass(unfold: easyunfold.unfold.UnfoldKSet, intensity_tol: float = 0.1, extrema_tol: float = 0.001, parabolic: bool = True, npoints: float = 3)
:canonical: easyunfold.effective_mass.EffectiveMass

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass.__init__
```

````{py:method} set_nocc(nocc)
:canonical: easyunfold.effective_mass.EffectiveMass.set_nocc

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass.set_nocc
```

````

````{py:property} kpoints
:canonical: easyunfold.effective_mass.EffectiveMass.kpoints

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass.kpoints
```

````

````{py:property} kpoints_labels
:canonical: easyunfold.effective_mass.EffectiveMass.kpoints_labels

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass.kpoints_labels
```

````

````{py:method} get_band_extrema(mode: str = 'cbm', extrema_tol: float = None, ispin=0)
:canonical: easyunfold.effective_mass.EffectiveMass.get_band_extrema

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass.get_band_extrema
```

````

````{py:method} _get_kpoint_distances()
:canonical: easyunfold.effective_mass.EffectiveMass._get_kpoint_distances

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass._get_kpoint_distances
```

````

````{py:method} _get_fitting_data(kidx: int, iband: int, direction=1, ispin=0, npoints=None)
:canonical: easyunfold.effective_mass.EffectiveMass._get_fitting_data

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass._get_fitting_data
```

````

````{py:method} get_effective_masses(npoints: typing.Union[float, None] = None, ispin=0, iks=None, iband=None, mode=None)
:canonical: easyunfold.effective_mass.EffectiveMass.get_effective_masses

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass.get_effective_masses
```

````

````{py:method} _get_effective_masses(mode: str = 'cbm', ispin: int = 0, npoints: typing.Union[None, int] = None, iks=None, iband=None)
:canonical: easyunfold.effective_mass.EffectiveMass._get_effective_masses

```{autodoc2-docstring} easyunfold.effective_mass.EffectiveMass._get_effective_masses
```

````

`````

````{py:function} locate_kpoint_segment(idxk: int, label_idx: list, label_names: list, direction: int)
:canonical: easyunfold.effective_mass.locate_kpoint_segment

```{autodoc2-docstring} easyunfold.effective_mass.locate_kpoint_segment
```
````
