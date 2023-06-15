"""
Module for obtaining effective mass
"""
from typing import Union

import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import curve_fit

from .unfold import UnfoldKSet
# pylint: disable=invalid-name

eV_to_hartree = physical_constants['electron volt-hartree relationship'][0]
bohr_to_m = physical_constants['Bohr radius'][0]
angstrom_to_bohr = bohr_to_m / 1e-10

TMP_DATA = {}


def fit_effective_mass(distances, energies, parabolic=True):
    """Fit the effective masses using either a parabolic or nonparabolic fit.

    Adapted from ``sumo``.

    Args:
        distances (:obj:`numpy.ndarray`): The x-distances between k-points in
            reciprocal Angstroms, normalised to the band extrema.
        energies (:obj:`numpy.ndarray`): The band eigenvalues normalised to the
            eigenvalue of the band extrema.
        parabolic (:obj:`bool`, optional): Use a parabolic fit of the band
            edges. If ``False`` then nonparabolic fitting will be attempted.
            Defaults to ``True``.

    Returns:
        float: The effective mass in units of electron rest mass, :math:`m_0`.
    """
    if parabolic:
        fit = np.polyfit(distances, energies, 2)
        c = 2 * fit[0]  # curvature therefore 2 * the exponent on the ^2 term

    else:
        # Use non parabolic description of the bands
        def f(x, alpha, d):
            top = np.sqrt(4 * alpha * d * x**2 + 1) - 1
            bot = 2 * alpha
            return top / bot

        # set boundaries for curve fitting: alpha > 1e-8
        # as alpha = 0 causes an error
        bounds = ((1e-8, -np.inf), (np.inf, np.inf))
        popt, _ = curve_fit(f, distances, energies, p0=[1.0, 1.0], bounds=bounds)  # pylint: disable=unbalanced-tuple-unpacking
        c = 2 * popt[1]

    # coefficient is currently in eV/Angstrom^2/h_bar^2
    # want it in atomic units so Hartree/bohr^2/h_bar^2
    eff_mass = (angstrom_to_bohr**2 / eV_to_hartree) / c
    return eff_mass


def fitted_band(x: np.ndarray, eff_mass: float) -> np.ndarray:
    """Return fitted effective mass curve"""
    c = (angstrom_to_bohr**2 / eV_to_hartree) / eff_mass
    x0 = x - x[0]
    return x0 + x[0], c / 2 * x0**2


def points_with_tol(array, value, tol=1e-4):
    """
    Return the indices and values of points in an array close to the value with a tolerance
    """
    idx = np.where(np.abs(array - value) < tol)[0]
    return idx, array[idx]


class EffectiveMass:
    """Calculate effective mass from unfolding data"""

    def __init__(
        self,
        unfold: UnfoldKSet,
        intensity_tol: float = 1e-1,
        extrema_tol: float = 1e-3,
        degeneracy_tol: float = 1e-2,
        parabolic: bool = True,
        npoints: float = 3,
    ):
        """
        Instantiate the object

        Args:
            unfold (UnfoldKSet): The ``UnfoldKSet`` object that holes unfolding data.
            intensity_tol (float): Intensity tolerance for detecting band edges.
            parabolic (bool): Perform parabolic fit or not. Defaults to True as non-parabolic fit is not working at the moment...
        """
        self.unfold: UnfoldKSet = unfold
        self.intensity_tol = intensity_tol
        self.extrema_detect_tol = extrema_tol
        self.degeneracy_tol = degeneracy_tol
        self.parabolic = parabolic
        self.nocc = None  # Number of occupied bands
        if npoints < 3:
            raise ValueError('At least three points are needed for fitting the effective mass!')

        self.npoints = npoints

    def set_nocc(self, nocc):
        self.nocc = nocc

    @property
    def kpoints(self):
        return self.unfold.kpts_pc

    @property
    def kpoints_labels(self):
        return self.unfold.kpoint_labels

    def get_band_extrema(self, mode: str = 'cbm', extrema_tol: float = None, degeneracy_tol: float = None, ispin=0):
        """
        Obtain the kpoint idx of band maximum, sub indices in th set and the band indices.

        The search takes two steps, first the kpoints at the band extrema is located by comparing the
        band energies with that recorded in supplied *cbm* and *vbm*, based on the `exgtrema_tol`.
        Afterwards, the band indices are selected at the these kpoints using `degeneracy_tol`.

        Returns:
            A tuple of extrema locations including a list of kpoint indices, sub-indices within
            the set and the band indices at each kpoint that is within the `tol` set.
        """
        if extrema_tol is None:
            extrema_tol = self.extrema_detect_tol
        if degeneracy_tol is None:
            degeneracy_tol = self.degeneracy_tol

        intensity_tol = self.intensity_tol

        if mode not in ['cbm', 'vbm']:
            raise ValueError(f'Unknown mode {mode}')
        cbm = self.unfold.calculated_quantities[mode]
        weights = self.unfold.calculated_quantities['spectral_weights_per_set']

        # Indices of the kpoint corresponding to the CBM
        k_indicies = []
        k_subset_indices = []
        cbm_indices = []
        for ik, wset in enumerate(weights):
            for isubset in range(wset.shape[1]):
                if np.any(np.abs(wset[ispin, isubset, :, 0] - cbm) < extrema_tol):
                    itmp, _ = points_with_tol(wset[ispin, isubset, :, 0], cbm, extrema_tol)
                    # Check if it has sufficient intensity
                    if np.max(wset[ispin, isubset, itmp, 1]) < intensity_tol:
                        continue
                    # Select this kpoints
                    k_indicies.append(ik)
                    # Select all band indices within the tolerance
                    k_subset_indices.append(isubset)
                    # Stop looking at other kpoints in the k subset if found
                    break

        # Go through each case
        for ik, iksub in zip(k_indicies, k_subset_indices):
            itmp, _ = points_with_tol(weights[ik][ispin, iksub, :, 0], cbm, degeneracy_tol)
            cbm_indices.append(itmp)

        return k_indicies, k_subset_indices, cbm_indices

    def _get_kpoint_distances(self):
        """
        Distances between the kpoints along the path in the reciprocal space.
        This does not take account of the breaking of the path
        NOTE: the reciprocal lattice vectors includes the 2pi factor, e.g. np.linalg.inv(L).T * 2 * np.pi
        """
        kpts = self.kpoints
        pc_latt = self.unfold.pc_latt

        kpts_path = kpts @ np.linalg.inv(pc_latt).T * np.pi * 2  # Kpoint path in the reciprocal space

        dists = np.cumsum(np.linalg.norm(np.diff(kpts_path, axis=0), axis=-1))
        dists = np.append([0], dists)
        return dists

    def _get_fitting_data(self, kidx: int, iband: int, direction=1, ispin=0, npoints=None):
        """
        Get fitting data for a specific combination of kpoint and band index
        """
        istart = kidx
        weights = self.unfold.calculated_quantities['spectral_weights_per_set']
        dists = self._get_kpoint_distances()
        kdists = []
        engs_effective = []

        npoints = self.get_npoints(npoints)
        for i in range(npoints):
            idx = istart + i * direction
            kdists.append(dists[idx])
            # Get the spectral weight array
            sw = weights[idx]
            intensities = sw[ispin, :, iband, 1]
            engs = sw[ispin, :, iband, 0]
            kw = self.unfold.expansion_results['weights'][idx]
            # Compute the effective energy weighted by intensity and kpoint weighting
            eng_effective = np.sum(engs * intensities * kw) / np.sum(intensities * kw)
            engs_effective.append(eng_effective)
        return kdists, engs_effective

    def get_npoints(self, override: Union[float, None] = None):
        """Get the number of points used for fitting"""
        if override is None:
            return self.npoints
        return override

    def get_effective_masses(self, npoints: Union[float, None] = None, ispin=0):
        """
        Workout the effective masses based on the unfolded band structure
        """
        outputs = {}
        for mode in ['cbm', 'vbm']:
            name = 'electrons' if mode == 'cbm' else 'holes'
            outputs[name] = self._get_effective_masses(mode, npoints=npoints, ispin=ispin)
        return outputs

    def _get_effective_masses(self, mode: str = 'cbm', ispin: int = 0, npoints: Union[None, int] = None):
        """
        Work out the effective masses based on the unfolded band structure for CBM or VBM
        """
        iks, _, iband = self.get_band_extrema(mode=mode)
        # Override occupations
        if self.nocc:
            iband = [self.nocc for _ in iband]

        npoints = self.get_npoints(npoints)
        results = []
        label_idx = [x[0] for x in self.kpoints_labels]
        label_names = [x[1] for x in self.kpoints_labels]

        for idxk, idxb in zip(iks, iband):
            for direction in (-1, 1):
                # Check if the direction makes sense
                if idxk + direction < 0 or idxk + direction >= len(self.kpoints):
                    continue
                # Check if the direction is broken
                if idxk + direction in label_idx:
                    continue
                # Get fitting data for each (degenerate) band at the extrema
                for band_id in idxb:
                    kdists, engs_effective = self._get_fitting_data(idxk, band_id, direction, ispin=ispin, npoints=npoints)

                    me = fit_effective_mass(kdists, engs_effective, parabolic=self.parabolic)

                    # If the identified edge is not in the list of high symmetry point, ignore it
                    # This mitigate the problem where the CBM can be duplicated....
                    ilabel, label_from, label_to = locate_kpoint_segment(idxk, label_idx, label_names, direction)
                    # Record the results
                    results.append({
                        'kpoint_index': idxk,
                        'direction': direction,
                        'effective_mass': me,
                        'band_index': band_id,
                        'kpoint_label_from': label_from,
                        'kpoint_from': self.kpoints[idxk],
                        'kpoint_label_to': label_to,
                        'kpoint_to': self.kpoints[label_idx[ilabel + direction]],
                        'type': 'electrons' if mode == 'cbm' else 'holes',
                        'raw_data': {
                            'kpoint_distances': kdists,
                            'effective_energies': engs_effective
                        }
                    })
        results.sort(key=lambda x: abs(abs(x['effective_mass'])))
        results.sort(key=lambda x: abs(x['kpoint_index']))
        return results


def locate_kpoint_segment(idxk: int, label_idx: list, label_names: list, direction: int):
    """Locate the labels and indices of the kpoints defining a segment"""
    if idxk not in label_idx:
        pairs = list(zip(label_idx, label_names))
        i = 0
        for i, (ikto, _) in enumerate(pairs):
            if ikto > idxk:
                break
        label_to = label_names[i]
        label_from = 'N/A'
        ilabel = label_idx.index(label_idx[i - 1])
    else:
        ilabel = label_idx.index(idxk)
        label_from = label_names[ilabel]
        label_to = label_names[ilabel + direction]
    return ilabel, label_from, label_to
