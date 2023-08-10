# {py:mod}`easyunfold.procar`

```{py:module} easyunfold.procar
```

```{autodoc2-docstring} easyunfold.procar
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Procar <easyunfold.procar.Procar>`
  - ```{autodoc2-docstring} easyunfold.procar.Procar
    :summary:
    ```
````

### API

`````{py:class} Procar(fobjs_or_paths=None, is_soc=False)
:canonical: easyunfold.procar.Procar

Bases: {py:obj}`monty.json.MSONable`

```{autodoc2-docstring} easyunfold.procar.Procar
```

```{rubric} Initialization
```

```{autodoc2-docstring} easyunfold.procar.Procar.__init__
```

````{py:method} _read(fobj, parsed_kpoints=None)
:canonical: easyunfold.procar.Procar._read

```{autodoc2-docstring} easyunfold.procar.Procar._read
```

````

````{py:method} _read_header_nion_proj_names(fobj)
:canonical: easyunfold.procar.Procar._read_header_nion_proj_names

```{autodoc2-docstring} easyunfold.procar.Procar._read_header_nion_proj_names
```

````

````{py:method} read(fobjs_or_paths)
:canonical: easyunfold.procar.Procar.read

```{autodoc2-docstring} easyunfold.procar.Procar.read
```

````

````{py:method} get_projection(atom_idx: typing.List[int], proj: typing.Union[typing.List[str], str], weight_by_k=False)
:canonical: easyunfold.procar.Procar.get_projection

```{autodoc2-docstring} easyunfold.procar.Procar.get_projection
```

````

````{py:method} as_dict() -> dict
:canonical: easyunfold.procar.Procar.as_dict

```{autodoc2-docstring} easyunfold.procar.Procar.as_dict
```

````

````{py:method} from_dict(d)
:canonical: easyunfold.procar.Procar.from_dict
:classmethod:

```{autodoc2-docstring} easyunfold.procar.Procar.from_dict
```

````

`````
