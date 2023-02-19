#


## EffectiveMass
```python 
EffectiveMass(
   unfold: UnfoldKSet, intensity_tol = 0.1, extrema_tol = 0.001, degeneracy_tol = 0.01,
   parabolic = True
)
```


---
Calculate effective mass from unfolding data


**Methods:**


### .set_nocc
```python
.set_nocc(
   nocc
)
```


### .kpoints
```python
.kpoints()
```


### .kpoints_labels
```python
.kpoints_labels()
```


### .get_band_extrema
```python
.get_band_extrema(
   mode: str = 'cbm', extrema_tol: float = None, degeneracy_tol: float = None, ispin = 0
)
```

---
Obtain the kpoint idx of band maximum, sub indices in th set and the band indices.

The search takes two steps, first the kpoints at the band extrema is located by comparing the
band energies with that recorded in supplied *cbm* and *vbm*, based on the `exgtrema_tol`.
Afterwards, the band indices are selected at the these kpoints using `degeneracy_tol`.


**Returns**

A tuple of extrema locations including a list of kpoint indices, sub-indices within
the set and the band indices at each kpoint that is within the `tol` set.

### .get_effective_masses
```python
.get_effective_masses(
   npoints = 3, ispin = 0
)
```

---
Workout the effective masses based on the unfolded band structure

----


### fit_effective_mass
```python
.fit_effective_mass(
   distances, energies, parabolic = True
)
```

---
Fit the effective masses using either a parabolic or nonparabolic fit.

Adapted from ``sumo``.


**Args**

* **distances** (:obj:`numpy.ndarray`) : The x-distances between k-points in
    reciprocal Angstroms, normalised to the band extrema.
* **energies** (:obj:`numpy.ndarray`) : The band eigenvalues normalised to the
    eigenvalue of the band extrema.
* **parabolic** (:obj:`bool`, optional) : Use a parabolic fit of the band
    edges. If ``False`` then nonparabolic fitting will be attempted.
    Defaults to ``True``.


**Returns**

* **float**  : The effective mass in units of electron rest mass, :math:`m_0`.


----


### points_with_tol
```python
.points_with_tol(
   array, value, tol = 0.0001
)
```

---
Return the indices and values of points in an array close to the value with a tolerance
