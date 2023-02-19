#


## UnfoldPlotter
```python 
UnfoldPlotter(
   unfold: UnfoldKSet
)
```


---
A collection plotting tools for unfolded band structures


**Methods:**


### .plot_spectral_function
```python
.plot_spectral_function(
   engs, sf, eref = None, figsize = (4, 3), ylim = (-3, 3), dpi = 150, vscale = 1.0,
   contour_plot = False, alpha = 1.0, save = False, ax = None, vmin = None, vmax = None,
   cmap = 'PuRd', show = False, title = None
)
```

---
Plot spectral function.


**Args**

* **engs** (numpy.ndarray) : The energies of the spectral functions.
* **sf** (np.ndarray) : An array of the spectral function.
* **eref** (float) : Reference energy to be used - this energy will be set as zero.
* **figsize**  : Size of the figure.
* **ylim**  : Plotting limit for the y-axis, with respect to `eref`.
* **dpi**  : DPI of the generated graph.
* **vscale**  : A scaling factor for the colour map.
* **contour_plot** (bool) : Whether to use contour plot instead of normal meshed color map.
* **alpha** (float) : Alpha for the color map.
* **save** (str) : Name of the file where the generated figure is saved.
* **ax**  : An existing plotting axis to be be used.
* **vmin** (float) : Minimum value for the color map.
* **vmax** (float) : Maximum value for the color map.
* **cmap** (str) : Name of the color map to be used.


**Returns**

The figure generated containing the spectral function.

### .plot_spectral_function_rgba
```python
.plot_spectral_function_rgba(
   engs, sf, eref = None, figsize = (4, 3), ylim = (-3, 3), dpi = 150, intensity = 1.0,
   save = False, ax = None, show = False, title = None, vmin = None, vmax = None
)
```

---
Plot spectral function.


**Args**

* **engs** (numpy.ndarray) : The energies of the spectral functions.
* **sf** (np.ndarray) : An array of the spectral function.
* **eref** (float) : Reference energy to be used - this energy will be set as zero.
* **figsize**  : Size of the figure.
* **ylim**  : Plotting limit for the y-axis, with respect to `eref`.
* **dpi**  : DPI of the generated graph.
* **vscale**  : A scaling factor for the colour map.
* **contour_plot** (bool) : Whether to use contour plot instead of normal meshed color map.
* **alpha** (float) : Alpha for the color map.
* **save** (str) : Name of the file where the generated figure is saved.
* **ax**  : An existing plotting axis to be be used.
* **vmin** (float) : Minimum value for the color map.
* **vmax** (float) : Maximum value for the color map.
* **cmap** (str) : Name of the color map to be used.


**Returns**

The figure generated containing the spectral function.

### .plot_effective_mass
```python
.plot_effective_mass(
   eff: EffectiveMass, engs, sf, eref = None, save = None, show = False,
   effective_mass_data = None, **kwargs
)
```

---
Plot the effective masses on top of the spectral function.


**Args**

* **eff**  : An `EffectiveMass` object used for plotting.
* **engs** (numpy.ndarray) : The energies of the spectral functions.
* **sf** (np.ndarray) : An array of the spectral function.
* **eref** (float) : Reference energy to be used - this energy will be set as zero.
* **kwargs**  : Other keyword arguments supplied to `plot_spectral_function`.


**Returns**

A figure with the data used for fitting effective mass plotted on top of the spectral function.

### .plot_spectral_weights
```python
.plot_spectral_weights(
   figsize = (4, 3), ylim = (-3, 3), dpi = 150, factor = 3.0, eref = None, color = 'C1',
   alpha = 0.5, save = None, ax = None, show = False, title = None
)
```

---
Plot the spectral weights.


**Note**

The reduction of symmetry means there can be multiple supercell kpoints for each
primitive cell kpoint. When using this scattering plot representation, the markers can
overlap each other leading to misrepresentations of the actual effective band structure.

---
However, this function is still useful when: 1. the symmetry splitting is turned off. 2.
direct visualization of the underlying spectral weight is needed. 3. Check the correctness
of effective mass extraction.


**Args**

* **eref** (float) : Reference energy to be used - this energy will be set as zero.
* **figsize**  : Size of the figure.
* **ylim**  : Plotting limit for the y-axis, with respect to `eref`.
* **dpi**  : DPI of the generated graph.
* **alpha** (float) : Alpha for the markers.
* **save** (str) : Name of the file where the generated figure is saved.
* **ax**  : Existing plotting axes to be be used (list if having two spin channels).
* **factor** (float) : Scaling factor for the marker size.
* **color** (str) : Color of the markers.


**Returns**

A Figure with the spectral weights plotted as a scatter plot.

### .plot_projected
```python
.plot_projected(
   procar, eref = None, gamma = False, npoints = 2000, sigma = 0.2, ncl = False,
   symm_average = True, figsize = (4, 3), ylim = (-3, 3), dpi = 150, vscale = 1.0,
   contour_plot = False, alpha = 1.0, save = False, ax = None, vmin = None, vmax = None,
   cmap = 'PuRd', show = False, title = None, atoms_idx = None, orbitals = None,
   intensity = 1.0, use_subplot = False, colors = None, colorspace = 'lab'
)
```

---
Plot projected sepctral function onto multiple subplots or a single plot with colormapping.

This simply computes the spectral function at each orbital/atoms sites and plot them onto
multiple subplots. The columns are for each orbital and the rows are for each spin channel.

----


### interpolate_colors
```python
.interpolate_colors(
   colors, weights, colorspace = 'lab', normalize = True
)
```

---
Interpolate colors at a number of points within a colorspace.


**Args**

* **colors** (str) : A list of colors specified in any way supported by matplotlib.
* **weights** (list) : A list of weights with the shape (n, N). Where the N values of
    the last axis give the amount of N colors supplied in `colors`.
* **colorspace** (str) : The colorspace in which to perform the interpolation. The
    allowed values are rgb, hsv, lab, luvlc, lablch, and xyz.


**Returns**

A list of colors, specified in the rgb format as a (n, 3) array.

----


### process_projection_options
```python
.process_projection_options(
   atoms_idx, orbitals
)
```

---
Process commandline type specifications
