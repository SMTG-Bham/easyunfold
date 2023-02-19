#


## UnfoldKSet
```python 
UnfoldKSet(
   M, kpts_pc, pc_latt, pc_opts, sc_opts, time_reversal = True, expand = True,
   metadata = None, expansion_results = None, calculated_quantities = None,
   kpoint_labels = None
)
```


---
Stores the information of the kpoints in the primitive cell, and what they unfolds to in the supercell


**Methods:**


### .is_calculated
```python
.is_calculated()
```

---
Show the status of the work

### .has_averaged_spectral_weights
```python
.has_averaged_spectral_weights()
```

---
Return True if the spectral weights stored is averaged

### .from_atoms
```python
.from_atoms(
   cls, M, kpts_pc, pc, sc, time_reversal = True, expand = True, symprec = 1e-05
)
```

---
Initialise from primitive cell and supercell atoms

### .expand_pc_kpoints
```python
.expand_pc_kpoints()
```

---
Comptue the pc kpoints to be unfolded into

### .nsymm_orig
```python
.nsymm_orig()
```

---
Number of symmetry operation in the original cell

### .nsymm_expand
```python
.nsymm_expand()
```

---
Number of symmetry operation in the original cell

### .nkpts_orig
```python
.nkpts_orig()
```

---
Total number of unexpanded kpoints

### .nkpts_expand
```python
.nkpts_expand()
```

---
Total number of expanded kpoints

### .generate_sc_kpoints
```python
.generate_sc_kpoints()
```

---
Generate the supercell kpoints to be calculated


**Returns**

A flat list of supercell kpoints in fractional coordinates
An indexing nested list to map expanded kpoints set to the supercell kpoints generated

### .write_sc_kpoints
```python
.write_sc_kpoints(
   file, nk_per_split = None, scf_kpoints_and_weights = None
)
```

---
Write the supercell kpoints

### .write_pc_kpoints
```python
.write_pc_kpoints(
   file, expanded = False
)
```

---
Write the primitive cell kpoints

### .load_procar
```python
.load_procar(
   procar: Union[str, List[str]], force = False
)
```

---
Read in PROCAR for band-based projection

### .procars
```python
.procars()
```

---
Loaded PROCARS

### .procar_kmaps
```python
.procar_kmaps()
```

---
Loaded PROCARS

### .get_band_weight_sets
```python
.get_band_weight_sets(
   atoms_idx, orbitals, procars = None
)
```

---
Get weights array sets for bands
Construct the weights of each band in same format of the kpoint set.
Each item is an numpy array of (nspins, nbands), containing the summed weights over
the passed atom indices and orbitals.

### .get_spectral_function
```python
.get_spectral_function(
   wavecar = None, npoints = 2000, sigma = 0.1, gamma = False, ncl = False,
   gamma_half = 'x', symm_average = True, atoms_idx = None, orbitals = None
)
```

---
Get the spectral function

### .get_spectral_weights
```python
.get_spectral_weights(
   wavecar = None, gamma = False, ncl = False, gamma_half = 'x', symm_average = True
)
```

---
Get the spectral function

### .as_dict
```python
.as_dict()
```

---
To a dictionary representation

### .get_kpoint_distances
```python
.get_kpoint_distances()
```

---
Distances between the kpoints along the path in the reciprocal space.
This does not take account of the breaking of the path.
NOTE: the reciprocal lattice vectors includes the 2pi factor, e.g. np.linalg.inv(L).T * 2 * np.pi

----


## Unfold
```python 
Unfold(
   M = None, wavecar = 'WAVECAR', gamma = False, lsorbit = False, gamma_half = 'x',
   verbose = False
)
```


---
Low lever interface for performing unfolding related calculations.
obtain the effective band structure (EBS).

REF:
"Extracting E versus k effective band structure from supercell
 calculations on alloys and impurities"
Phys. Rev. B 85, 085201 (2012)


**Methods:**


### .get_vbm_cbm
```python
.get_vbm_cbm(
   thresh = 1e-08
)
```

---
Locate the VBM from the WAVECAR

### .get_ovlap_G
```python
.get_ovlap_G(
   ikpt = 1, epsilon = 1e-05
)
```

---
Get subset of the reciprocal space vectors of the supercell,
specifically the ones that match the reciprocal space vectors of the
primitive cell.

### .find_K_index
```python
.find_K_index(
   K0
)
```

---
Find the index of K0.

### .k2K_map
```python
.k2K_map(
   kpath
)
```

---
Find the map from primitive-cell k-points to supercell k-points.

### .spectral_weight_k
```python
.spectral_weight_k(
   k0, whichspin = 1
)
```

---
Spectral weight for a given k:

P_{Km}(k) = \sum_n |<Km | kn>|^2

---
which is equivalent to

    P_{Km}(k) = \sum_{G} |C_{Km}(G + k - K)|^2

where {G} is a subset of the reciprocal space vectors of the supercell.

### .spectral_weight
```python
.spectral_weight(
   kpoints
)
```

---
Calculate the spectral weight for a list of kpoints in the primitive BZ.

### .spectral_function
```python
.spectral_function(
   nedos = 4000, sigma = 0.02
)
```

---
Generate the spectral function

A(k_i, E) = \sum_m P_{Km}(k_i)\Delta(E - Em)

---
Where the \Delta function can be approximated by Lorentzian or Gaussian
function.

----


### get_symmetry_dataset
```python
.get_symmetry_dataset(
   atoms, symprec = 1e-05
)
```

---
Get the symmetry dataset using spglib

----


### find_K_from_k
```python
.find_K_from_k(
   k: np.ndarray, M: np.ndarray
)
```

---
Get the K vector of the supercell onto which the k vector of the primitive
cell folds. The unfolding vector G, which satisfy the following equation,
is also returned.

k = K + G

---
where G is a reciprocal space vector of supercell

----


### rotate_kpt
```python
.rotate_kpt(
   k: np.ndarray, opt: np.ndarray
)
```

---
Apply rotation to a kpoint based on the rotations of the crystals (in the real space)

NOTE: The rotation matrix should be the one that act on fractional coordinates, e.g. from spglib.

----


### expand_K_by_symmetry
```python
.expand_K_by_symmetry(
   kpt, opts_pc, opts_sc, time_reversal = True
)
```

---
Expend the sampling of the PC kpoints due to symmetry breaking of the SC

----


### LorentzSmearing
```python
.LorentzSmearing(
   x, x0, sigma = 0.02
)
```

---
Simulate the Delta function by a Lorentzian shape function

\Delta(x) = \lim_{\sigma\to 0}  Lorentzian

----


### GaussianSmearing
```python
.GaussianSmearing(
   x, x0, sigma = 0.02
)
```

---
Simulate the Delta function by a Lorentzian shape function

\Delta(x) = \lim_{\sigma\to 0} Gaussian

----


### removeDuplicateKpoints
```python
.removeDuplicateKpoints(
   kpoints, return_map = False, decimals = 6
)
```

---
remove duplicate kpoints in the list.

----


### write_kpoints
```python
.write_kpoints(
   kpoints: Union[np.ndarray, list], outpath = 'KPOINTS', comment = '', weights = None
)
```

---
save to VASP KPOINTS file

----


### read_kpoints
```python
.read_kpoints(
   path = 'KPOINTS'
)
```

---
Read kpoints from a KPOINTS file containing reciprocal space coordinates (fractional)

Returns the kpoints, the comment and the labels at each kpoint

----


### read_kpoints_line
```python
.read_kpoints_line(
   content, density = 20
)
```

---
Read kpoints in the line mode

Resolve to explicit kpoints

----


### make_kpath
```python
.make_kpath(
   kbound, nseg = 40
)
```

---
Return a list of kpoints defining the path between the given kpoints.

----


### EBS_scatter
```python
.EBS_scatter(
   kpts, cell, spectral_weight, atomic_weights = None, atomic_colors = None,
   eref = 0.0, nseg = None, save = 'ebs_s.png', kpath_label = None, factor = 20,
   figsize = (3.0, 4.0), ylim = (-3, 3), show = True, ax = None, color = 'b'
)
```

---
plot the effective band structure with scatter, the size of the scatter
indicates the spectral weight.
The plotting function utilizes Matplotlib package.

inputs:
kpts: the kpoints vectors in fractional coordinates.
cell: the primitive cell basis
spectral_weight: self-explanatory

----


### EBS_cmaps
```python
.EBS_cmaps(
   kpts, cell, E0, spectral_function, eref = 0.0, nseg = None, kpath_label = None,
   explicit_labels = None, save = None, figsize = (3.0, 4.0), ylim = (-3, 3), show = True,
   contour_plot = False, ax = None, vscale = 1.0, title = None, vmax = None, vmin = None,
   alpha = 1.0, cmap = 'jet'
)
```

---
plot the effective band structure with colormaps.  The plotting function
utilizes Matplotlib package.


**Args**

* **kpts**  : the kpoints vectors in fractional coordinates.
* **cell**  : the primitive cell basis
* **E0**  : The energies corresponds to each element of the spectral function
* **spectral_function**  : The spectral function array in the shape of (nspin, nk, neng)
* **eref**  : Refernce point for zero energy
* **kpath_label**  : Label of the high symmetry kpoints along the pathway
* **nseg**  : Number of points in each segment of the kpoint pathway
* **explicit_labels**  : A list of tuplies containing tuples of `(index, label)` to explicitly set kpoint labels.
* **save**  : Name of the file the plot to be saved to.
* **figsize**  : Size of hte figure
* **ylim**  : Limit for the y axis. The limit is applied *after* substracting the refence energy.
* **show**  : To show the plot interactively or not.
* **contour_plot**  : Plot in the contour mode
* **ax**  : Existing axis(axes) to plot onto
* **cmap**  : Colour mapping for the density/contour plot
* **title**  : Title to be used
* **vscale**  : Scale factor for color coding


----


### clean_latex_string
```python
.clean_latex_string(
   label
)
```

---
Clean up latex labels and convert if necessary

----


### spectral_function_from_weight_sets
```python
.spectral_function_from_weight_sets(
   spectral_weight_sets, kweight_sets, nedos = 4000, sigma = 0.02, emin = None,
   emax = None, band_weight_sets = None
)
```

---
Generate the spectral function

A(k_i, E) = \sum_m P_{Km}(k_i)\Delta(E - Em)

---
Where the \Delta function can be approximated by Lorentzian or Gaussian
function.


**Args**

* **band_weight_sets** (np.ndarray) : Additional weighting for each band, used for generating
  projection onto atomic orbitals.

----


### spectral_weight_multiple_source
```python
.spectral_weight_multiple_source(
   kpoints, unfold_objs, transform_matrix
)
```

---
Calculate the spectral weight for a list of kpoints in the primitive BZ
from a list of WAVECARs.

----


### concatenate_scf_kpoints
```python
.concatenate_scf_kpoints(
   scf_kpts, scf_weights, kpoints
)
```

---
Concatenate SCF kpoints (from IBZKPT) with zero-weighted kpoints

----


### create_white_colormap
```python
.create_white_colormap(
   color: Union[str, tuple, list]
)
```

---
Create colormap from white to certain colour.


**Args**

* **color** (str, tuple, list) : HEX color string or tuple/list of RGB color.


----


### create_white_colormap_from_existing
```python
.create_white_colormap_from_existing(
   name: str
)
```

---
Create a white-based color map from an existing one.

----


### parse_atoms_idx
```python
.parse_atoms_idx(
   atoms_idx
)
```

---
Expanding syntex like `1-2` (inclusive)

For example, `1,2,3,4-6` will be expanded as `1,2,3,4,5,6`.
