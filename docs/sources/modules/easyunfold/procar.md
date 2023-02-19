#


## Procar
```python 
Procar(
   fobj_or_path = None, is_soc = False
)
```


---
Reader for PROCAR file


**Methods:**


### .get_projection
```python
.get_projection(
   atom_idx: List[int], proj: Union[List[str], str], weight_by_k = False
)
```

---
Get project for specific atoms and specific projectors


**Args**

* **atom_idx** (list) : A list of index of the atoms to be selected
* **proj** (list) : A list of the projector names to be selected
* **weight_by_k**  : Apply k weighting or not.


**Returns**

The project summed over the selected atoms and the projectors
