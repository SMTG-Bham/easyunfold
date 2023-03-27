# {py:mod}`easyunfold.wavecar`

```{py:module} easyunfold.wavecar
```

```{autodoc2-docstring} easyunfold.wavecar
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Wavecar <easyunfold.wavecar.Wavecar>`
  - ```{autodoc2-docstring} easyunfold.wavecar.Wavecar
    :summary:
    ```
````

### API

`````{py:class} Wavecar(fnm='WAVECAR', lsorbit=False, lgamma=False, gamma_half='x')
:canonical: easyunfold.wavecar.Wavecar

```{autodoc2-docstring} easyunfold.wavecar.Wavecar
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.__init__
```

````{py:method} close()
:canonical: easyunfold.wavecar.Wavecar.close

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.close
```

````

````{py:method} is_soc()
:canonical: easyunfold.wavecar.Wavecar.is_soc

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.is_soc
```

````

````{py:method} is_gamma()
:canonical: easyunfold.wavecar.Wavecar.is_gamma

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.is_gamma
```

````

````{py:method} read_wf_header()
:canonical: easyunfold.wavecar.Wavecar.read_wf_header

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.read_wf_header
```

````

````{py:method} set_wf_prec()
:canonical: easyunfold.wavecar.Wavecar.set_wf_prec

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.set_wf_prec
```

````

````{py:method} read_bands() -> tuple
:canonical: easyunfold.wavecar.Wavecar.read_bands

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.read_bands
```

````

````{py:method} get_gvectors(ikpt=1, force_gamma=False, check_consistency=True)
:canonical: easyunfold.wavecar.Wavecar.get_gvectors

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.get_gvectors
```

````

````{py:method} read_band_coeffs(ispin=1, ikpt=1, iband=1, norm=False)
:canonical: easyunfold.wavecar.Wavecar.read_band_coeffs

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.read_band_coeffs
```

````

````{py:method} locate_rec(ispin=1, ikpt=1, iband=1)
:canonical: easyunfold.wavecar.Wavecar.locate_rec

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.locate_rec
```

````

````{py:method} check_index(ispin, ikpt, iband)
:canonical: easyunfold.wavecar.Wavecar.check_index

```{autodoc2-docstring} easyunfold.wavecar.Wavecar.check_index
```

````

`````
