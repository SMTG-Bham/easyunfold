# {py:mod}`easyunfold.utils`

```{py:module} easyunfold.utils
```

```{autodoc2-docstring} easyunfold.utils
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`write_kpoints <easyunfold.utils.write_kpoints>`
  - ```{autodoc2-docstring} easyunfold.utils.write_kpoints
    :summary:
    ```
* - {py:obj}`write_kpoints_castep <easyunfold.utils.write_kpoints_castep>`
  - ```{autodoc2-docstring} easyunfold.utils.write_kpoints_castep
    :summary:
    ```
* - {py:obj}`write_kpoints_vasp <easyunfold.utils.write_kpoints_vasp>`
  - ```{autodoc2-docstring} easyunfold.utils.write_kpoints_vasp
    :summary:
    ```
* - {py:obj}`read_kpoints_line_vasp <easyunfold.utils.read_kpoints_line_vasp>`
  - ```{autodoc2-docstring} easyunfold.utils.read_kpoints_line_vasp
    :summary:
    ```
* - {py:obj}`read_kpoints <easyunfold.utils.read_kpoints>`
  - ```{autodoc2-docstring} easyunfold.utils.read_kpoints
    :summary:
    ```
* - {py:obj}`read_kpoints_vasp <easyunfold.utils.read_kpoints_vasp>`
  - ```{autodoc2-docstring} easyunfold.utils.read_kpoints_vasp
    :summary:
    ```
* - {py:obj}`read_kpoints_castep <easyunfold.utils.read_kpoints_castep>`
  - ```{autodoc2-docstring} easyunfold.utils.read_kpoints_castep
    :summary:
    ```
* - {py:obj}`wrap_kpoints <easyunfold.utils.wrap_kpoints>`
  - ```{autodoc2-docstring} easyunfold.utils.wrap_kpoints
    :summary:
    ```
* - {py:obj}`find_unique <easyunfold.utils.find_unique>`
  - ```{autodoc2-docstring} easyunfold.utils.find_unique
    :summary:
    ```
* - {py:obj}`reduce_kpoints <easyunfold.utils.reduce_kpoints>`
  - ```{autodoc2-docstring} easyunfold.utils.reduce_kpoints
    :summary:
    ```
* - {py:obj}`kpoints_equal <easyunfold.utils.kpoints_equal>`
  - ```{autodoc2-docstring} easyunfold.utils.kpoints_equal
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`RE_COMMENT <easyunfold.utils.RE_COMMENT>`
  - ```{autodoc2-docstring} easyunfold.utils.RE_COMMENT
    :summary:
    ```
````

### API

````{py:data} RE_COMMENT
:canonical: easyunfold.utils.RE_COMMENT
:value: >
   None

```{autodoc2-docstring} easyunfold.utils.RE_COMMENT
```

````

````{py:function} write_kpoints(kpoints: typing.Union[numpy.ndarray, list], outpath, *args, code='vasp', **kwargs)
:canonical: easyunfold.utils.write_kpoints

```{autodoc2-docstring} easyunfold.utils.write_kpoints
```
````

````{py:function} write_kpoints_castep(kpoints: typing.Union[numpy.ndarray, list], dest, source=None, tag='spectral', weights=None, **kwargs)
:canonical: easyunfold.utils.write_kpoints_castep

```{autodoc2-docstring} easyunfold.utils.write_kpoints_castep
```
````

````{py:function} write_kpoints_vasp(kpoints: typing.Union[numpy.ndarray, list], outpath: str = 'KPOINTS', comment: str = '', weights: typing.Union[None, typing.List[float]] = None, **kwargs)
:canonical: easyunfold.utils.write_kpoints_vasp

```{autodoc2-docstring} easyunfold.utils.write_kpoints_vasp
```
````

````{py:function} read_kpoints_line_vasp(content, density=20)
:canonical: easyunfold.utils.read_kpoints_line_vasp

```{autodoc2-docstring} easyunfold.utils.read_kpoints_line_vasp
```
````

````{py:function} read_kpoints(path='KPOINTS', code='vasp', **kwargs)
:canonical: easyunfold.utils.read_kpoints

```{autodoc2-docstring} easyunfold.utils.read_kpoints
```
````

````{py:function} read_kpoints_vasp(path='KPOINTS')
:canonical: easyunfold.utils.read_kpoints_vasp

```{autodoc2-docstring} easyunfold.utils.read_kpoints_vasp
```
````

````{py:function} read_kpoints_castep(path, tag='spectral')
:canonical: easyunfold.utils.read_kpoints_castep

```{autodoc2-docstring} easyunfold.utils.read_kpoints_castep
```
````

````{py:function} wrap_kpoints(kpoints: typing.Union[list, numpy.ndarray])
:canonical: easyunfold.utils.wrap_kpoints

```{autodoc2-docstring} easyunfold.utils.wrap_kpoints
```
````

````{py:function} find_unique(seq: numpy.ndarray, func=None)
:canonical: easyunfold.utils.find_unique

```{autodoc2-docstring} easyunfold.utils.find_unique
```
````

````{py:function} reduce_kpoints(kpoints: typing.Union[list, numpy.ndarray], time_reversal=True, rounding_digits=10)
:canonical: easyunfold.utils.reduce_kpoints

```{autodoc2-docstring} easyunfold.utils.reduce_kpoints
```
````

````{py:function} kpoints_equal(k1, k2, time_reversal=False, atol=1e-05)
:canonical: easyunfold.utils.kpoints_equal

```{autodoc2-docstring} easyunfold.utils.kpoints_equal
```
````
