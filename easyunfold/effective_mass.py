"""
Module for obtaining effective mass
"""
import numpy as np
from scipy.constants import physical_constants
from scipy.optimize import curve_fit

from .unfold import UnfoldKSet
# pylint: disable=invalid-name

eV_to_hartree = physical_constants['electron volt-hartree relationship'][0]
bohr_to_m = physical_constants['Bohr radius'][0]
angstrom_to_bohr = bohr_to_m / 1e-10


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
        popt, _ = curve_fit(f, distances, energies, p0=[1.0, 1.0], bounds=bounds)
        c = 2 * popt[1]

    # coefficient is currently in eV/Angstrom^2/h_bar^2
    # want it in atomic units so Hartree/bohr^2/h_bar^2
    eff_mass = (angstrom_to_bohr**2 / eV_to_hartree) / c
    return eff_mass


def points_with_tol(array, value, tol=1e-4):
    """
    Return the indices and values of points in an array close to the value with a tolerance
    """
    idx = np.where(np.abs(array - value) < tol)[0]
    return idx, array[idx]


class EffectiveMass:
    """Calculate effective mass from unfolding data"""

    def __init__(self, unfold: UnfoldKSet, intensity_tol=1e-1, edge_detect_tol=1e-2, parabolic=True):
        """
        Instantiate the object

        Args:
            unfold (UnfoldKSet): The ``UnfoldKSet`` object that holes unfolding data.
            intensity_tol (float): Intensity tolerance for detecting band edges.
            parabolic (bool): Perform parabolic fit or not. Defaults to True as non-parabolic fit is not working at the moment...
        """
        self.unfold: UnfoldKSet = unfold
        self.intensity_tol = intensity_tol
        self.edge_detect_tol = edge_detect_tol
        self.parabolic = parabolic
        self.nocc = None  # Number of occupied bands

    def set_nocc(self, nocc):
        self.nocc = nocc

    @property
    def kpoints(self):
        return self.unfold.kpts_pc

    @property
    def kpoints_labels(self):
        return self.unfold.kpoint_labels

    def get_band_extrema(self, mode: str = 'cbm', tol: float = None):
        """
        Obtain the kpoint idx of band maximum

        Returns:
            A tuple of extrema locations including a list of kpoint indices, sub-indices within
            the set and the band indices.
        """
        if tol is None:
            tol = self.edge_detect_tol
        intensity_tol = self.intensity_tol

        if mode not in ['cbm', 'vbm']:
            raise ValueError(f'Unknown mode {mode}')
        cbm = self.unfold.calculated_quantities[mode]
        weights = self.unfold.calculated_quantities['spectral_weights_per_set']

        # Indices of the kpoint corresponding to the CBM
        k_indicies = []
        k_subset_indices = []
        cbm_indices = []
        for idx, wset in enumerate(weights):
            tmp_spin = wset[0]
            for isubset, tmp in enumerate(tmp_spin):
                if np.any(np.abs(tmp[:, 0] - cbm) < tol):
                    itmp, _ = points_with_tol(tmp[:, 0], cbm, tol)
                    # Check if it has sufficient intensity
                    if tmp[min(itmp), 1] < intensity_tol:
                        continue
                    k_indicies.append(idx)

                    if mode == 'cbm':
                        cbm_indices.append(min(itmp))
                    else:
                        cbm_indices.append(max(itmp))
                    k_subset_indices.append(isubset)
                    break

        if np.any(min(cbm_indices) != cbm_indices):
            print(f'WARNING: band max detection failure - detected band indices: {cbm_indices}')
        return k_indicies, k_subset_indices, cbm_indices

    def _get_kpoint_distances(self):
        """
        Distances between the kpoints along the path in the reciprocal space.
        This does not take account of the breaking of the path
        NOTE: the reciprocal lattice vectors includes the 2pi factor, e.g. np.linalg.inv(L) * 2 * np.pi
        """
        kpts = self.kpoints
        pc_latt = self.unfold.pc_latt

        kpts_path = kpts @ np.linalg.inv(pc_latt) * np.pi * 2  # Kpoint path in the reciprocal space

        dists = np.cumsum(np.linalg.norm(np.diff(kpts_path, axis=0), axis=-1))
        dists = np.append([0], dists)
        return dists

    def _get_fitting_data(self, kidx: int, iband: int, direction=1, ispin=0, npoints=3):
        """
        Get fitting data for a specific combination of kpoint and band index
        """
        istart = kidx
        weights = self.unfold.calculated_quantities['spectral_weights_per_set']
        dists = self._get_kpoint_distances()
        kdists = []
        engs_effective = []
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

    def get_effective_masses(self, npoints=3, ispin=0):
        """
        Workout the effective masses based on the unfolded band structure
        """
        outputs = {}
        for mode in ['cbm', 'vbm']:
            name = 'electrons' if mode == 'cbm' else 'holes'
            outputs[name] = self._get_effective_masses(mode, npoints=npoints, ispin=ispin)
        return outputs

    def _get_effective_masses(self, mode='cbm', ispin=0, npoints=3):
        """
        Workout the effective masses based on the unfolded band structure for CBM or VBM
        """
        iks, _, iband = self.get_band_extrema(mode=mode)
        # Override occupations
        if self.nocc:
            iband = [self.nocc for _ in iband]

        results = []
        label_idx = [x[0] for x in self.kpoints_labels]
        label_names = [x[1] for x in self.kpoints_labels]

        for idxk in iks:
            for direction in (-1, 1):
                # Check if the direction makes sense
                if idxk + direction < 0 or idxk + direction >= len(self.kpoints):
                    continue
                # Check if the direction is broken
                if idxk + direction in label_idx:
                    continue
                # Get fitting data
                kdists, engs_effective = self._get_fitting_data(idxk, iband[0], direction, ispin=ispin, npoints=npoints)
                me = fit_effective_mass(kdists, engs_effective, parabolic=self.parabolic)

                # If the identified edge is not in the list of high symmetry point, ignore it
                # This mitigate the problem where the CBM can be duplicated....
                if idxk not in label_idx:
                    continue

                ilabel = label_idx.index(idxk)
                label_from = label_names[ilabel]
                label_to = label_names[ilabel + direction]

                # Record the results
                results.append({
                    'kpoint_index': idxk,
                    'direction': direction,
                    'effective_mass': me,
                    'kpoint_label_from': label_from,
                    'kpoint_from': self.kpoints[idxk],
                    'kpoint_label_to': label_to,
                    'kpoint_to': self.kpoints[label_idx[ilabel + direction]],
                    'type': 'electrons' if mode == 'cbm' else 'holes'
                })
        results.sort(key=lambda x: abs(x['effective_mass']))
        return results
