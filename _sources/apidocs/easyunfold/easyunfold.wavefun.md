# {py:mod}`easyunfold.wavefun`

```{py:module} easyunfold.wavefun
```

```{autodoc2-docstring} easyunfold.wavefun
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`WaveFunction <easyunfold.wavefun.WaveFunction>`
  - ```{autodoc2-docstring} easyunfold.wavefun.WaveFunction
    :summary:
    ```
* - {py:obj}`VaspWaveFunction <easyunfold.wavefun.VaspWaveFunction>`
  - ```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction
    :summary:
    ```
* - {py:obj}`CastepWaveFunction <easyunfold.wavefun.CastepWaveFunction>`
  - ```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction
    :summary:
    ```
````

### API

`````{py:class} WaveFunction(wfc)
:canonical: easyunfold.wavefun.WaveFunction

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.__init__
```

````{py:property} kpoints
:canonical: easyunfold.wavefun.WaveFunction.kpoints
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.kpoints
```

````

````{py:property} nkpts
:canonical: easyunfold.wavefun.WaveFunction.nkpts

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.nkpts
```

````

````{py:property} nspins
:canonical: easyunfold.wavefun.WaveFunction.nspins
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.nspins
```

````

````{py:property} mesh_size
:canonical: easyunfold.wavefun.WaveFunction.mesh_size
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.mesh_size
```

````

````{py:property} bands
:canonical: easyunfold.wavefun.WaveFunction.bands
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.bands
```

````

````{py:property} nbands
:canonical: easyunfold.wavefun.WaveFunction.nbands

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.nbands
```

````

````{py:property} occupancies
:canonical: easyunfold.wavefun.WaveFunction.occupancies
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.occupancies
```

````

````{py:method} get_gvectors(ik)
:canonical: easyunfold.wavefun.WaveFunction.get_gvectors
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.get_gvectors
```

````

````{py:method} get_band_coeffs(ispin, ik, ib, norm=True)
:canonical: easyunfold.wavefun.WaveFunction.get_band_coeffs
:abstractmethod:

```{autodoc2-docstring} easyunfold.wavefun.WaveFunction.get_band_coeffs
```

````

`````

`````{py:class} VaspWaveFunction(wfc: easyunfold.wavecar.Wavecar)
:canonical: easyunfold.wavefun.VaspWaveFunction

Bases: {py:obj}`easyunfold.wavefun.WaveFunction`

```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction.__init__
```

````{py:property} kpoints
:canonical: easyunfold.wavefun.VaspWaveFunction.kpoints

````

````{py:property} occupancies
:canonical: easyunfold.wavefun.VaspWaveFunction.occupancies

````

````{py:property} nspins
:canonical: easyunfold.wavefun.VaspWaveFunction.nspins

````

````{py:property} mesh_size
:canonical: easyunfold.wavefun.VaspWaveFunction.mesh_size

```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction.mesh_size
```

````

````{py:property} bands
:canonical: easyunfold.wavefun.VaspWaveFunction.bands

```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction.bands
```

````

````{py:method} get_gvectors(ik)
:canonical: easyunfold.wavefun.VaspWaveFunction.get_gvectors

```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction.get_gvectors
```

````

````{py:method} get_band_coeffs(ispin, ik, ib, norm=True)
:canonical: easyunfold.wavefun.VaspWaveFunction.get_band_coeffs

```{autodoc2-docstring} easyunfold.wavefun.VaspWaveFunction.get_band_coeffs
```

````

`````

`````{py:class} CastepWaveFunction(wfc: castepxbin.wave.WaveFunction)
:canonical: easyunfold.wavefun.CastepWaveFunction

Bases: {py:obj}`easyunfold.wavefun.WaveFunction`

```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction.__init__
```

````{py:method} from_file(fname)
:canonical: easyunfold.wavefun.CastepWaveFunction.from_file
:classmethod:

```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction.from_file
```

````

````{py:property} kpoints
:canonical: easyunfold.wavefun.CastepWaveFunction.kpoints

````

````{py:property} occupancies
:canonical: easyunfold.wavefun.CastepWaveFunction.occupancies

````

````{py:property} mesh_size
:canonical: easyunfold.wavefun.CastepWaveFunction.mesh_size

```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction.mesh_size
```

````

````{py:property} bands
:canonical: easyunfold.wavefun.CastepWaveFunction.bands

````

````{py:property} nspins
:canonical: easyunfold.wavefun.CastepWaveFunction.nspins

````

````{py:method} get_gvectors(ik)
:canonical: easyunfold.wavefun.CastepWaveFunction.get_gvectors

```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction.get_gvectors
```

````

````{py:method} get_band_coeffs(ispin, ik, ib, norm=True)
:canonical: easyunfold.wavefun.CastepWaveFunction.get_band_coeffs

```{autodoc2-docstring} easyunfold.wavefun.CastepWaveFunction.get_band_coeffs
```

````

`````
