#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The main module for unfolding workflow and algorithm
"""
# pylint: disable=invalid-name,protected-access,too-many-locals

############################################################
from typing import Union, List, Tuple
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from monty.json import MSONable
from tqdm import tqdm
import spglib

from easyunfold import __version__
from .pyvaspwfc.vaspwfc import vaspwfc
from .pyvaspwfc.procar import procar

############################################################


def get_symmetry_dataset(atoms, symprec=1e-5):
    """Get the symmetry dataset using spglib"""
    return spglib.get_symmetry_dataset((atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers()), symprec=symprec)


def find_K_from_k(k: np.ndarray, M: np.ndarray):
    """
    Get the K vector of the supercell onto which the k vector of the primitive
    cell folds. The unfolding vector G, which satisfy the following equation,
    is also returned.

        k = K + G

    where G is a reciprocal space vector of supercell
    """

    M = np.array(M)
    Kc = np.dot(k, M.T)
    G = np.array(np.round(Kc), dtype=int)
    KG = Kc - np.round(Kc)

    return KG, G


def rotate_kpt(k: np.ndarray, opt: np.ndarray):
    """
    Apply rotation to a kpoint based on the rotations of the crystals (in the real space)

    NOTE: The rotation matrix should be the one that act on fractional coordinates, e.g. from spglib.
    """
    rot = k @ opt
    return rot - np.rint(rot)


def expand_K_by_symmetry(kpt, opts_pc, opts_sc, time_reversal=True):
    """
    Expend the sampling of the PC kpoints due to symmetry breaking of the SC
    """
    kpt_orig = np.asarray(kpt)

    # Find distinct images of the kpoints in the PC
    pc_distinct = [kpt_orig]
    for opt in opts_pc:
        k_equiv = rotate_kpt(kpt_orig, opt)
        k_equiv_neg = k_equiv * -1
        found = False
        for k_ in pc_distinct:
            if np.allclose(k_equiv, k_):
                found = True
            if time_reversal and np.allclose(k_equiv_neg, k_):
                found = True
        if not found:
            pc_distinct.append(k_equiv)
    weights = np.ones(len(pc_distinct))
    # Now each for uniqueness in the supercell
    for i, kpt_ in enumerate(pc_distinct):
        if weights[i] == 0:
            continue

        for opt in opts_sc:
            k_equiv = rotate_kpt(kpt_, opt)
            k_equiv_neg = k_equiv * -1
            for j, k_ in enumerate(pc_distinct):
                # Skip if it is the same point, or the other point has been taken
                if i == j or weights[j] == 0:
                    continue
                if np.allclose(k_equiv, k_) or (time_reversal and np.allclose(k_equiv_neg, k_)):
                    weights[i] += weights[j]
                    weights[j] = 0

    out_points = []
    out_weights = []
    for i, kpt_ in enumerate(pc_distinct):
        if weights[i] != 0:
            out_points.append(kpt_)
            out_weights.append(weights[i])
    assert sum(out_weights) == len(pc_distinct)
    # The first kpoint of the set should always be the original kpoint
    assert out_points[0] is kpt_orig
    out_weights = np.array(out_weights) / sum(out_weights)
    return out_points, out_weights


class UnfoldKSet(MSONable):
    """Stores the information of the kpoints in the primitive cell, and what they unfolds to in the supercell"""

    _VERSION = '0.1.0'

    def __init__(self,
                 M,
                 kpts_pc,
                 pc_latt,
                 pc_opts,
                 sc_opts,
                 time_reversal=True,
                 expand=True,
                 metadata=None,
                 expansion_results=None,
                 calculated_quantities=None,
                 kpoint_labels=None):
        """
        Args:
            kpts_pc: A list of kpoints in the PC
            pc_lattice: A 3x3 matrix of row lattice vectors of the primitive cell
            pc_opts: Symmetry operations of the primitive cell
            sc_opts: Symmetry operations of the supercell
            expand: Expand the kpoint to take account of broken symmetry or not
        """
        # Basic properties - needed to recreate the object
        self.kpts_pc = kpts_pc
        self.pc_latt = pc_latt
        self.pc_opts = pc_opts
        self.sc_opts = sc_opts
        self.expand = expand
        self.M = M
        self.expansion_results = expansion_results
        self.time_reversal = time_reversal
        self.calculated_quantities = {} if not calculated_quantities else calculated_quantities
        self.kpoint_labels = kpoint_labels
        if metadata is None:
            metadata = {}
        self.metadata = metadata
        self.metadata['program_version'] = __version__

        # Transient properties
        self.reduced_sckpts = None
        self.reduced_sckpts_map = None
        if self.expansion_results is None:
            self.expand_pc_kpoints()

    @property
    def is_calculated(self):
        """Show the status of the work"""
        return bool(self.calculated_quantities)

    @property
    def has_averaged_spectral_weights(self):
        """Return True if the spectral weights stored is averaged"""
        return self.calculated_quantities.get('spectral_weights_is_averaged', False)

    @classmethod
    def from_atoms(cls, M, kpts_pc, pc, sc, time_reversal=True, expand=True, symprec=1e-5):
        """Initialise from primitive cell and supercell atoms"""
        pc_symm_data = get_symmetry_dataset(pc, symprec=symprec)
        sc_symm_data = get_symmetry_dataset(sc, symprec=symprec)
        return cls(
            M=M,
            kpts_pc=kpts_pc,
            pc_latt=np.asarray(pc.cell),
            pc_opts=pc_symm_data['rotations'],
            sc_opts=sc_symm_data['rotations'],
            time_reversal=time_reversal,
            expand=expand,
            metadata={
                'symmetry_dataset_pc': pc_symm_data,
                'symmetry_dataset_sc': sc_symm_data,
            },
        )

    def expand_pc_kpoints(self):
        """Comptue the pc kpoints to be unfolded into"""
        expended_k = []
        expended_weights = []
        for kpt in self.kpts_pc:
            # Record the expanded kpoints and corresponding weights
            if self.expand:
                kset, weights = expand_K_by_symmetry(kpt, self.pc_opts, self.sc_opts, time_reversal=self.time_reversal)
            else:
                # Just take the original point and set the weight to be unity
                kset = [kpt]
                weights = np.array([1.0])
            expended_k.append(kset)
            expended_weights.append(weights)
        # Save as the attribute
        self.expansion_results = {'kpoints': expended_k, 'weights': expended_weights}

    def __repr__(self) -> str:
        return f'<UnfoldKSet with {self.nkpts_expand}/{self.nkpts_orig} kpoints based on {self.nsymm_expand}/{self.nsymm_orig} symm ops'

    @property
    def nsymm_orig(self):
        """Number of symmetry operation in the original cell"""
        return self.pc_opts.shape[0]

    @property
    def nsymm_expand(self):
        """Number of symmetry operation in the original cell"""
        return self.sc_opts.shape[0]

    @property
    def nkpts_orig(self):
        """Total number of unexpanded kpoints"""
        return len(self.expansion_results['kpoints'])

    @property
    def nkpts_expand(self):
        """Total number of expanded kpoints"""
        return sum(map(len, self.expansion_results['kpoints']))

    def generate_sc_kpoints(self):
        """
        Generate the supercell kpoints to be calculated

        Returns:
            A flat list of supercell kpoints in fractional coordinates
            An indexing nested list to map expanded kpoints set to the supercell kpoints generated
        """

        assert bool(self.expansion_results)
        expended_sc = []
        all_sc = []
        # Find the SC kpoint for each PC kpoint
        for kset in self.expansion_results['kpoints']:
            this_k = []
            for kpt in kset:
                sc_k, _ = find_K_from_k(kpt, self.M)
                this_k.append(sc_k)
            expended_sc.append(this_k)
            all_sc.extend(this_k)
        # We now have bunch of supercell kpoints for each set of expanded kpoints
        # Try to find duplicated SC kpoints
        reduced_sckpts, sc_kpts_map = removeDuplicateKpoints(all_sc, return_map=True)
        sc_kpts_map = list(sc_kpts_map)

        # Mapping between the pckpts to the reduced_sckpts
        reduced_sc_map = []
        for sc_set in expended_sc:
            map_indx = []
            for _ in sc_set:
                map_indx.append(sc_kpts_map.pop(0))
            reduced_sc_map.append(map_indx)

        self.expansion_results['reduced_sckpts'] = reduced_sckpts
        self.expansion_results['reduced_sckpts_map'] = reduced_sc_map
        # A nested list that stores the indices of the sc kpts in the reduced_sckpts list

    def write_sc_kpoints(self, file, nk_per_split=None, scf_kpoints_and_weights=None):
        """Write the supercell kpoints"""
        if self.expansion_results.get('reduced_sckpts') is None:
            self.generate_sc_kpoints()
        kpoints = np.asarray(self.expansion_results['reduced_sckpts'])
        weights = None
        if nk_per_split is None:
            if scf_kpoints_and_weights:
                # Prepend with SCF kpoints
                kpoints, weights = concatenate_scf_kpoints(scf_kpoints_and_weights[0], scf_kpoints_and_weights[1], kpoints)

            write_kpoints(kpoints, file, comment='supercell kpoints', weights=weights)
        else:
            splits = [kpoints[i:i + nk_per_split] for i in range(0, kpoints.shape[0], nk_per_split)]
            for i_spilt, kpt in enumerate(splits):
                if scf_kpoints_and_weights:
                    kpt, weights = concatenate_scf_kpoints(scf_kpoints_and_weights[0], scf_kpoints_and_weights[1], kpt)
                write_kpoints(kpt, f'{file}_{i_spilt + 1:03d}', f'supercell kpoints split {i_spilt + 1}', weights)

    def write_pc_kpoints(self, file, expanded=False):
        """Write the primitive cell kpoints"""
        if expanded:
            all_pc = []
            for tmp in self.expansion_results['kpoints']:
                all_pc.extend(tmp)
            write_kpoints(all_pc, file, comment='expanded primitive cell kpoints')
        else:
            write_kpoints(self.kpts_pc, file, comment='primitive cell kpoints')

    def _read_weights(self, wavecar: Union[str, List[str]], gamma: bool, ncl: bool, gamma_half: str):
        """
        Read the weights from the wavecar for all kpoints

        Returns the averaged weights and the original weights per set of kpoints
        """
        weights_per_set = []
        averaged_weights = []
        if not isinstance(wavecar, (list, tuple)):
            wavecar = [wavecar]
        # Load the WAVECAR unfold objects
        unfold_objs = [Unfold(self.M, name, gamma=gamma, lsorbit=ncl, gamma_half=gamma_half) for name in wavecar]
        # Record the VBM and the CBM values
        varray = np.array([obj.get_vbm_cbm() for obj in unfold_objs])
        self.calculated_quantities['vbm'] = float(varray[:, 0].max())
        self.calculated_quantities['cbm'] = float(varray[:, 1].min())
        # Read the spectral weights for each set of expanded kpoints
        for kset, weights in tqdm(zip(self.expansion_results['kpoints'], self.expansion_results['weights']),
                                  desc='kpt',
                                  total=len(self.expansion_results['kpoints'])):
            # Raw spectral weights
            sw = spectral_weight_multiple_source(kset, unfold_objs, self.M)
            weights_per_set.append(sw.copy())
            # Take weighted average
            for ik, w in enumerate(weights):
                sw[:, ik, :, :] *= w
            averaged_weights.append(sw.sum(axis=1))

        # Recreate the full weights array
        averaged_weights = np.stack(averaged_weights, axis=1)
        self.calculated_quantities['spectral_weights_per_set'] = weights_per_set
        self.calculated_quantities['version'] = self._VERSION

        return averaged_weights, weights_per_set

    def _get_spectral_weights(self,
                              wavecar,
                              npoints=2000,
                              sigma=0.01,
                              emin=None,
                              emax=None,
                              gamma=False,
                              ncl=False,
                              gamma_half='x',
                              also_spectral_function=False,
                              symm_average=True):
        """
        Fetch spectral weights from a wavecar and compute spectral function is requested
        """
        # If WAVECAR is given - reload from the data
        if wavecar:
            self._read_weights(wavecar, gamma=gamma, ncl=ncl, gamma_half=gamma_half)
        else:
            if not self.is_calculated:
                raise RuntimeWarning('The spectral weights need to be calculated first - please pass the WAVECAR file(s).')

        # Use existing results
        if symm_average:
            sws = self.calculated_quantities['spectral_weights_per_set']
            kweight_sets = self.expansion_results['weights']
        else:
            # No averaging - we just return the first item of each set, which is the weight of the original set
            sws = [item[:, :1, :, :] for item in self.calculated_quantities['spectral_weights_per_set']]
            # Recreate the full weights array
            kweight_sets = [[1.0] for i in range(len(sws))]

        if also_spectral_function:
            e0, spectral_function = spectral_function_from_weight_sets(sws, kweight_sets, nedos=npoints, sigma=sigma, emin=emin, emax=emax)
            self.calculated_quantities['e0'] = e0
            self.calculated_quantities['spectral_function'] = spectral_function
            return sws, e0, spectral_function
        return sws

    def get_spectral_function(self, wavecar=None, npoints=2000, sigma=0.02, gamma=False,
                              ncl=False, gamma_half='x', symm_average=True):
        """Get the spectral function"""
        _, e0, spectral_function = self._get_spectral_weights(wavecar,
                                                              npoints=npoints,
                                                              sigma=sigma,
                                                              gamma=gamma,
                                                              ncl=ncl,
                                                              gamma_half=gamma_half,
                                                              also_spectral_function=True,
                                                              symm_average=symm_average)
        return e0, spectral_function

    def get_spectral_weights(self, wavecar=None, npoints=2000, sigma=0.02, gamma=False, ncl=False,
                             gamma_half='x', symm_average=True):
        """Get the spectral function"""
        return self._get_spectral_weights(wavecar=wavecar,
                                          npoints=npoints,
					                      sigma=sigma,
                                          gamma=gamma,
                                          ncl=ncl,
                                          gamma_half=gamma_half,
                                          also_spectral_function=False,
                                          symm_average=symm_average)

    def as_dict(self):
        """To a dictionary representation"""
        output = {'@module': self.__class__.__module__, '@class': self.__class__.__name__, '@version': __version__}
        for key in [
                'M', 'kpts_pc', 'pc_latt', 'pc_opts', 'sc_opts', 'expansion_results', 'time_reversal', 'calculated_quantities',
                'kpoint_labels', 'expand', 'metadata'
        ]:
            output[key] = getattr(self, key)
        return output

    def get_reduced_atomic_weights(self, procar_path, natoms, ncl=False):
        p = procar(procar_path, lsoc=ncl) # parse PROCAR with pyvaspwfc
        kmap = self.expansion_results['reduced_sckpts_map']

        i = 0
        atomic_weights = []
        for iatom in tqdm(natoms, desc='Atom', total=len(natoms)):
            atomic_weights += [p.get_pw(f"{i:.0f}:{i+iatom:.0f}")]
            i += iatom
        atomic_weights = np.array(atomic_weights)

        reduced_atomic_weights = []
        for i in range(atomic_weights.shape[0]):  # per atom
            atom_awht = []
            for spin in range(atomic_weights.shape[1]):  # per spin
                spin_awht = []
                for pckpt in kmap:  # corresponding SC kpoints per PC kpoint
                    awhts = [atomic_weights[i, spin, sckpt, :] for sckpt in pckpt]
                    mean_awhts = np.array([np.mean(weight_tuple) for weight_tuple in zip(*awhts)])
                    # getting average atomic weights over corresponding SC kpoints for each PC
                    # kpoint, per band
                    spin_awht.append(mean_awhts)
                atom_awht.append(spin_awht)
            reduced_atomic_weights.append(atom_awht)

        self.calculated_quantities['atomic_weights'] = atomic_weights
        self.calculated_quantities['reduced_atomic_weights'] = reduced_atomic_weights

def LorentzSmearing(x, x0, sigma=0.02):
    r"""
    Simulate the Delta function by a Lorentzian shape function

        \Delta(x) = \lim_{\sigma\to 0}  Lorentzian
    """

    return 1. / np.pi * sigma**2 / ((x - x0)**2 + sigma**2)


def GaussianSmearing(x, x0, sigma=0.02):
    r"""
    Simulate the Delta function by a Lorentzian shape function

        \Delta(x) = \lim_{\sigma\to 0} Gaussian
    """

    return 1. / (np.sqrt(2 * np.pi) * sigma) * np.exp(-(x - x0)**2 / (2 * sigma**2))


def removeDuplicateKpoints(kpoints, return_map=False, decimals=6):
    """
    remove duplicate kpoints in the list.
    """
    kpoints = np.asarray(kpoints)
    _, kid, inv_kid = np.unique(
        np.round(kpoints, decimals),
        axis=0,
        return_index=True,
        return_inverse=True,
    )
    reducedK = kpoints[kid]

    if return_map:
        return reducedK, inv_kid
    return reducedK


def write_kpoints(kpoints: Union[np.ndarray, list], outpath='KPOINTS', comment='', weights=None):
    """
    save to VASP KPOINTS file
    """
    kpoints = np.asarray(kpoints)
    nkpts = kpoints.shape[0]

    with open(outpath, 'w') as vaspkpt:
        vaspkpt.write(comment + '\n')
        vaspkpt.write('%d\n' % nkpts)
        vaspkpt.write('Rec\n')
        for ik in range(nkpts):
            if weights is not None:
                line = '  %12.8f %12.8f %12.8f %12.8f\n' % (kpoints[ik, 0], kpoints[ik, 1], kpoints[ik, 2], weights[ik])
            else:
                line = '  %12.8f %12.8f %12.8f 1.0\n' % (kpoints[ik, 0], kpoints[ik, 1], kpoints[ik, 2])
            vaspkpt.write(line)


def read_kpoints(path='KPOINTS'):
    """
    Read kpoints from a KPOINTS file containing reciprocal space coordinates (fractional)

    Returns the kpoints, the comment and the labels at each kpoint
    """
    content = Path(path).read_text().split('\n')
    comment = content[0]
    nkpts = int(content[1])
    if content[2].lower().startswith('lin'):
        return read_kpoints_line(content)
    assert content[2].lower().startswith('rec'), 'Only Reciprocal space KPOINT file is supported'
    kpts = []
    labels = []
    weights = []
    ik = 0
    for line in content[3:]:
        tokens = line.split()
        this_kpt = [float(value) for value in tokens[:3]]
        weights.append(float(tokens[3]))

        if len(tokens) >= 5:
            labels.append([ik, tokens[4]])

        kpts.append(this_kpt)
        ik += 1
        if ik == nkpts:
            break
    return kpts, comment, labels, weights


def read_kpoints_line(content, density=20):
    """
    Read kpoints in the line mode

    Resolve to explicit kpoints
    """
    comment = content[0]
    density = int(content[1]) if content[1] else density

    assert content[2].lower().startswith('lin'), 'Only Line mode KPOINT file is supported'
    assert content[3].lower().startswith('rec'), 'Only Reciprocal coorindates are supported!'

    segs = []
    labels = []
    for line in content[4:]:
        tokens = line.split()
        if not tokens:
            continue
        point = [float(x) for x in tokens[:3]]
        labels.append(tokens[-1])
        segs.append(point)
    # Process the segments
    kpoints = []
    labels_loc = []
    for i in range(int(len(segs) / 2)):
        k1 = segs[i * 2]
        k2 = segs[i * 2 + 1]
        # Check for duplicate end point
        if kpoints and kpoints[-1] == k1:
            kpoints.pop()
            labels_loc.pop()
        this_seg = np.linspace(k1, k2, density)
        labels_loc.append((len(kpoints), labels[i]))
        kpoints.extend(this_seg.tolist())
        labels_loc.append((len(kpoints) - 1, labels[i + 1]))
    return kpoints, comment, labels_loc


def make_kpath(kbound, nseg=40):
    """
    Return a list of kpoints defining the path between the given kpoints.
    """
    kbound = np.array(kbound, dtype=float)
    kdist = np.diff(kbound, axis=0)

    kpath = [kbound[ii] + kdist[ii] / nseg * nk for ii in range(len(kdist)) for nk in range(nseg)]
    kpath.append(kbound[-1])
    return kpath


def EBS_scatter(kpts,
                cell,
                spectral_weight,
                atomic_weights=None,
                eref=0.0,
                nseg=None,
                save='ebs_s.png',
                kpath_label=None,
                factor=20,
                figsize=(3.0, 4.0),
                ylim=(-3, 3),
                show=True,
                ax=None,
                color='b'):
    """
    plot the effective band structure with scatter, the size of the scatter
    indicates the spectral weight.
    The plotting function utilizes Matplotlib package.

    inputs:
        kpts: the kpoints vectors in fractional coordinates.
        cell: the primitive cell basis
        spectral_weight: self-explanatory
    """

    atomic_colors = [] if atomic_colors is None else atomic_colors
    kpath_label = [] if kpath_label is None else kpath_label

    nspin = spectral_weight.shape[0]
    kpt_c = np.dot(kpts, np.linalg.inv(cell).T)
    kdist = np.r_[0, np.cumsum(np.linalg.norm(np.diff(kpt_c, axis=0), axis=1))]
    nb = spectral_weight.shape[2]
    x0 = np.tile(kdist, (nb, 1)).T

    if atomic_weights is not None:
        atomic_weights = np.asarray(atomic_weights)
        assert atomic_weights.shape[1:] == spectral_weight.shape[:-1]

        if not atomic_colors:
            atomic_colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']

    if ax is None:
        fig = plt.figure()
        if nspin == 1:
            axes = [plt.subplot(111)]
            fig.set_size_inches(figsize)
        else:
            axes = [plt.subplot(121), plt.subplot(122)]
            fig.set_size_inches((figsize[0] * 2, figsize[1]))
    else:
        if not isinstance(ax, list):
            axes = [ax]
        else:
            axes = ax
        fig = axes[0].figure

    for ispin in range(nspin):
        ax = axes[ispin]
        if atomic_weights is not None:
            for iatom in range(atomic_weights.shape[0]):
                ax.scatter(x0,
                           spectral_weight[ispin, :, :, 0] - eref,
                           s=spectral_weight[ispin, :, :, 1] * factor * atomic_weights[iatom][ispin, :, :],
                           lw=0.0,
                           color=atomic_colors[iatom])
        else:
            ax.scatter(x0, spectral_weight[ispin, :, :, 0] - eref, s=spectral_weight[ispin, :, :, 1] * factor, lw=0.0, color=color)

        ax.set_xlim(0, kdist.max())
        ax.set_ylim(*ylim)
        ax.set_ylabel('Energy [eV]', labelpad=5)

        if nseg:
            for kb in kdist[::nseg]:
                ax.axvline(x=kb, lw=0.5, color='k', ls=':', alpha=0.8)

            if kpath_label:
                ax.set_xticks(kdist[::nseg])
                kname = [x.upper() for x in kpath_label]
                for ii, _ in enumerate(kname):
                    if kname[ii] == 'G':
                        kname[ii] = r'$\mathrm{\mathsf{\Gamma}}$'
                    else:
                        kname[ii] = r'$\mathrm{\mathsf{%s}}$' % kname[ii]
                ax.set_xticklabels(kname)

    fig.tight_layout(pad=0.2)
    fig.savefig(save, dpi=360)
    if show:
        fig.show()


def EBS_cmaps(kpts,
              cell,
              E0,
              spectral_function,
              atomic_weights=None,
              atomic_energies=None,
              eref=0.0,
              nseg=None,
              kpath_label=None,
              explicit_labels=None,
              save=None,
              figsize=(3.0, 4.0),
              ylim=(-3, 3),
              show=True,
              contour_plot=False,
              ax=None,
              vscale=1.0,
              cmap='jet',): # TODO: add option to set atomic colours (other than rgb) and
    # interpolate between
    """
    plot the effective band structure with colormaps.  The plotting function
    utilizes Matplotlib package.

    Args:
        kpts: the kpoints vectors in fractional coordinates.
        cell: the primitive cell basis
        e0: The energies corresponds to each element of the spectral function
        spectral_function: The spectral function array in the shape of (nspin, nk, neng)
        atomic_weights: The atomic weights in the shape of (n_atoms, nspin, nk, neng)
        atomic_energies: The energies matching the atomic weights in the shape of (nspin, nk, neng)
        eref: Reference point for zero energy
        kpath_label: Label of the high symmetry kpoints along the pathway
        nseg: Number of points in each segment of the kpoint pathway
        explicit_labels: A list of tuples containing tuples of `(index, label)` to explicitly set kpoint labels.
        save: Name of the file the plot to be saved to.
        figsize: SIze of hte figure
        ylim: Limit for the y axis. The limit is applied *after* subtracting the reference energy.
        show: To show the plot interactively or not.
        contour_plot: Plot in the contour mode
        ax: Existing axis(axes) to plot onto
        vscale: Scale factor for color coding
        cmap: Colour mapping for the density/contour plot
    """

    kpath_label = [] if not kpath_label else kpath_label
    nspin = spectral_function.shape[0]
    kpt_c = np.dot(kpts, np.linalg.inv(cell).T)
    kdist = np.r_[0., np.cumsum(np.linalg.norm(np.diff(kpt_c, axis=0), axis=1))]
    xmin, xmax = kdist.min(), kdist.max()

    if ax is None:
        fig = plt.figure()
        if nspin == 1:
            axes = [plt.subplot(111)]
            fig.set_size_inches(figsize)
        else:
            axes = [plt.subplot(121), plt.subplot(122)]
            fig.set_size_inches((figsize[0] * 2, figsize[1]))
    else:
        if not isinstance(ax, list):
            axes = [ax]
        else:
            axes = ax
        fig = axes[0].figure

    X, Y = np.meshgrid(kdist, E0 - eref)

    # Calculate the min and max values within the field of view, scaled by the factor
    mask = (E0 - eref < ylim[1]) & (E0 - eref > ylim[0])  # non-critical missing "eref" in original
    vmax = spectral_function[:, mask, :].max()
    vmin = spectral_function[:, mask, :].min()
    vmax = (vmax - vmin) * vscale + vmin

    if atomic_weights is not None and atomic_energies is not None: # prepare atomic_weights RGB
        # array for plotting
        natoms = atomic_weights.shape[0]
        if natoms not in [2, 3]: # TODO: Add option to select 3 species to project from a set of
            # >3 in the structure
            raise ValueError("Atomic projections only supported for 2 or 3 atoms")
        weights = np.array([atomic_weights[iatom, :, :, :] for iatom in range(natoms)])  # natom,
        # nspin, nk, neng
        if natoms == 2: # just use red and blue
            weights = np.insert(weights, 1, np.zeros(weights[0].shape), axis=0)
        reshaped_weights = np.transpose(weights, (1, 2, 3, 0))  # nspin, nk, neng, natom

        rgbs = np.zeros((nspin, len(kdist), len(E0), 4))  # nspin, nk, ne, 4 (for RGBA given by
        # atomic_weights + spectral density
        for ispin in range(nspin):
            for ik, kdistance in enumerate(kdist):
                for ie, energy in enumerate(E0):
                    idx = (np.abs(atomic_energies[ispin, ik,
                                  :] - energy)).argmin()  # get matching band index for energy and k
                    rgbs[ispin, ik, ie, :] = np.append(reshaped_weights[ispin, ik, idx, :],  # RGB
                                                       1.5 * ((spectral_function[
                                                                   ispin, ie, ik] - vmin) / vmax))
                    # Alpha determined by spectral density

    for ispin in range(nspin):
        ax = axes[ispin]
        if contour_plot:
            ax.contourf(X, Y, spectral_function[ispin], cmap=cmap, vmax=vmax, vmin=vmin)
        elif atomic_weights is not None and atomic_energies is not None:
            ax.imshow(rgbs[ispin, :, mask, :], aspect='auto', origin='lower',  # nk, neng, 4
                          extent=[kdist[0], kdist[-1], ylim[0], ylim[-1]])
        else:
            ax.pcolormesh(X, Y, spectral_function[ispin], cmap=cmap, shading='auto', vmax=vmax, vmin=vmin)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(*ylim)
        ax.set_ylabel('Energy (eV)', labelpad=5)

        if nseg:
            for kb in kdist[::nseg]:
                ax.axvline(x=kb, lw=0.5, color='k', ls=':', alpha=0.8)

            if kpath_label:
                tick_pos = list(kdist[::nseg])
                ax.set_xticks(tick_pos)
                kname = [x.upper() for x in kpath_label]
                for ii, _ in enumerate(kname):
                    kname[ii] = clean_latex_string(kname[ii])
                ax.set_xticklabels(kname)
        elif explicit_labels:
            tick_locs = []
            tick_labels = []
            for index, label in explicit_labels:
                ax.axvline(x=kdist[index], lw=0.5, color='k', ls=':', alpha=0.8)
                tick_locs.append(kdist[index])
                tick_labels.append(label)
            ticks = []
            tick_labels = []
            for index, label in explicit_labels:
                ticks.append(kdist[index])
                tick_labels.append(clean_latex_string(label))
            ax.set_xticks(tick_locs)
            ax.set_xticklabels(tick_labels)

    fig.tight_layout(pad=0.2)
    if save:
        fig.savefig(save, dpi=300)
    if show:
        fig.show()
    return fig


def clean_latex_string(label):
    """Clean up latex labels and convert if necessary"""
    if label == 'G':
        label = r'$\mathrm{\mathsf{\Gamma}}$'
    elif label.startswith('\\'):  ## This is a latex formatted label already
        label = f'$\\mathrm{{\\mathsf{{{label}}}}}$'
    else:
        label = r'$\mathrm{\mathsf{%s}}$' % label
    return label


def spectral_function_from_weight_sets(spectral_weight_sets, kweight_sets, nedos=4000, sigma=0.02, emin=None, emax=None):
    r"""
    Generate the spectral function

        A(k_i, E) = \sum_m P_{Km}(k_i)\Delta(E - Em)

    Where the \Delta function can be approximated by Lorentzian or Gaussian
    function.
    """

    nk = len(spectral_weight_sets)
    ns = spectral_weight_sets[0].shape[0]
    spectral_function = np.zeros((ns, nedos, nk), dtype=float)

    emin = spectral_weight_sets[0][:, :, :, 0].min() if emin is None else emin
    emax = spectral_weight_sets[0][:, :, :, 0].max() if emax is None else emax
    e0 = np.linspace(emin - 5 * sigma, emax + 5 * sigma, nedos)

    for ispin in range(ns):
        for ii in range(nk):
            for jj in range(spectral_weight_sets[ii].shape[1]):
                kweight = kweight_sets[ii][jj]
                E_Km = spectral_weight_sets[ii][ispin, jj, :, 0]
                P_Km = spectral_weight_sets[ii][ispin, jj, :, 1]
                # Take weighted average spectral functions
                spectral_function[ispin, :,
                                  ii] += np.sum(LorentzSmearing(e0[:, np.newaxis], E_Km[np.newaxis, :], sigma=sigma) * P_Km[np.newaxis, :],
                                                axis=1) * kweight
    return e0, spectral_function


############################################################


class Unfold():
    """
    Low lever interface for performing unfolding related calculations.
    obtain the effective band structure (EBS).

    REF:
    "Extracting E versus k effective band structure from supercell
     calculations on alloys and impurities"
    Phys. Rev. B 85, 085201 (2012)
    """

    def __init__(self, M=None, wavecar='WAVECAR', gamma=False, lsorbit=False, gamma_half='x', verbose=False):
        """
        Initialization.

        M is the transformation matrix between supercell and primitive cell:

            M = np.dot(A, np.linalg.inv(a))

        In real space, the basis vectors of Supercell (A) and those of the
        primitive cell (a) satisfy:

            A = np.dot(M, a);      a = np.dot(np.linalg.inv(M), A)

        Whereas in reciprocal space

            b = np.dot(M.T, B);    B = np.dot(np.linalg.inv(M).T, b)

        wavecar is the location of VASP WAVECAR file that contains the
        wavefunction information of a supercell calculation.
        """

        # Whether the WAVECAR is a gamma-only version
        self._lgam = gamma
        self._lsoc = lsorbit

        self.M = np.array(M, dtype=float)
        assert self.M.shape == (3, 3), 'Shape of the tranformation matrix must be (3,3)'

        self.wfc = vaspwfc(wavecar, lsorbit=self._lsoc, lgamma=self._lgam, gamma_half=gamma_half)
        # all the K-point vectors
        self.kvecs = self.wfc._kvecs
        # all the KS energies in shape (ns, nk, nb)
        self.bands = self.wfc._bands

        # G-vectors within the cutoff sphere, let's just do it once for all.
        # self.allGvecs = np.array([self.wfc.gvectors(ikpt=kpt+1)
        #                           for kpt in range(self.wfc._nkpts)], dtype=int)

        # spectral weight for all the kpoints
        self.SW = None

        self.verbose = verbose

    def get_vbm_cbm(self, thresh=1e-8) -> Tuple[float, float]:
        """Locate the VBM from the WAVECAR"""
        occ = self.wfc._occs

        occupied = np.abs(occ) > thresh
        vbm = float(self.bands[occupied].max())
        cbm = float(self.bands[~occupied].min())
        return vbm, cbm

    def get_ovlap_G(self, ikpt=1, epsilon=1E-5):
        """
        Get subset of the reciprocal space vectors of the supercell,
        specifically the ones that match the reciprocal space vectors of the
        primitive cell.
        """

        assert 1 <= ikpt <= self.wfc._nkpts, 'Invalid K-point index!'

        # Reciprocal space vectors of the supercell in fractional unit
        Gvecs = self.wfc.gvectors(ikpt=ikpt)
        # Gvecs = self.allGvecs[ikpt - 1]

        if self._lgam:
            nplw = Gvecs.shape[0]
            tmp = np.zeros((nplw * 2 - 1, 3), dtype=int)
            # the gvectors of Gamma version only contains half the number of a
            # normal version.
            tmp[:nplw, ...] = Gvecs
            tmp[nplw:, ...] = -Gvecs[1:, ...]  # G' = -G

            Gvecs = tmp

        # Shape of Gvecs: (nplws, 3)
        # iGvecs = np.arange(Gvecs.shape[0], dtype=int)

        # Reciprocal space vectors of the primitive cell
        gvecs = np.dot(Gvecs, np.linalg.inv(self.M).T)
        # Deviation from the perfect sites
        gd = gvecs - np.round(gvecs)
        # match = np.linalg.norm(gd, axis=1) < epsilon
        match = np.alltrue(np.abs(gd) < epsilon, axis=1)

        # return Gvecs[match], iGvecs[match]
        return Gvecs[match], Gvecs

    def find_K_index(self, K0):
        """
        Find the index of K0.
        """

        for ii in range(self.wfc._nkpts):
            if np.alltrue(np.abs(self.wfc._kvecs[ii] - K0) < 1E-5):
                return ii + 1
        raise ValueError('Cannot find the corresponding K-points in WAVECAR!')

    def k2K_map(self, kpath):
        """
        Find the map from primitive-cell k-points to supercell k-points.
        """

        return [self.find_K_index(find_K_from_k(k, self.M)[0]) - 1 for k in kpath]

    def spectral_weight_k(self, k0, whichspin=1):
        r"""
        Spectral weight for a given k:

            P_{Km}(k) = \sum_n |<Km | kn>|^2

        which is equivalent to

            P_{Km}(k) = \sum_{G} |C_{Km}(G + k - K)|^2

        where {G} is a subset of the reciprocal space vectors of the supercell.
        """
        if self.verbose:
            print('Processing k-point %8.4f %8.4f %8.4f' % (k0[0], k0[1], k0[2]))

        # find the K0 onto which k0 folds
        # k0 = G0 + K0
        K0, G0 = find_K_from_k(k0, self.M)
        # find index of K0
        ikpt = self.find_K_index(K0)

        # get the overlap G-vectors
        Gvalid, Gall = self.get_ovlap_G(ikpt=ikpt)
        # Gnew = Gvalid + k0 - K0
        Goffset = Gvalid + G0[np.newaxis, :]

        # Index of the Gvalid in 3D grid
        GallIndex = Gall % self.wfc._ngrid[np.newaxis, :]
        GoffsetIndex = Goffset % self.wfc._ngrid[np.newaxis, :]

        # 3d grid for planewave coefficients
        wfc_k_3D = np.zeros(self.wfc._ngrid, dtype=np.complex128)

        # if self._lsoc:
        #     # the weights and corresponding energies
        #     P_Km = np.zeros((2, self.wfc._nbands), dtype=float)
        #     E_Km = np.zeros((2, self.wfc._nbands), dtype=float)
        # else:
        #     # the weights and corresponding energies
        #     P_Km = np.zeros(self.wfc._nbands, dtype=float)
        #     E_Km = np.zeros(self.wfc._nbands, dtype=float)

        # the weights and corresponding energies
        P_Km = np.zeros(self.wfc._nbands, dtype=float)
        E_Km = np.zeros(self.wfc._nbands, dtype=float)

        for nb in range(self.wfc._nbands):
            # initialize the array to zero, which is unnecessary since the
            # GallIndex is the same for the same K-point
            # wfc_k_3D[:,:,:] = 0.0

            if self._lsoc:
                # pad the coefficients to 3D grid
                band_coeff = self.wfc.readBandCoeff(ispin=whichspin, ikpt=ikpt, iband=nb + 1, norm=True)
                nplw = band_coeff.shape[0] // 2
                band_spinor_coeff = [band_coeff[:nplw], band_coeff[nplw:]]

                # energy
                E_Km[nb] = self.bands[whichspin - 1, ikpt - 1, nb]
                for Ispinor in range(2):
                    # band = band_spinor_coeff[Ispinor]
                    # band /= np.linalg.norm(band)
                    wfc_k_3D[GallIndex[:, 0], GallIndex[:, 1], GallIndex[:, 2]] = band_spinor_coeff[Ispinor]

                    # spectral weight
                    P_Km[nb] += np.linalg.norm(wfc_k_3D[GoffsetIndex[:, 0], GoffsetIndex[:, 1], GoffsetIndex[:, 2]])**2
            else:
                # pad the coefficients to 3D grid
                band_coeff = self.wfc.readBandCoeff(ispin=whichspin, ikpt=ikpt, iband=nb + 1, norm=True)
                if self._lgam:
                    nplw = band_coeff.size
                    tmp = np.zeros((nplw * 2 - 1), dtype=band_coeff.dtype)
                    # for Gamma version, the coefficients corresponding to G \ne 0
                    # is multiplied by a factor of sqrt(2)
                    band_coeff[1:] /= np.sqrt(2.)
                    tmp[:nplw] = band_coeff
                    tmp[nplw:] = band_coeff[1:].conj()
                    band_coeff = tmp

                wfc_k_3D[GallIndex[:, 0], GallIndex[:, 1], GallIndex[:, 2]] = band_coeff
                # energy
                E_Km[nb] = self.bands[whichspin - 1, ikpt - 1, nb]
                # spectral weight
                P_Km[nb] = np.linalg.norm(wfc_k_3D[GoffsetIndex[:, 0], GoffsetIndex[:, 1], GoffsetIndex[:, 2]])**2

        return np.array((E_Km, P_Km), dtype=float).T

    def spectral_weight(self, kpoints):
        """
        Calculate the spectral weight for a list of kpoints in the primitive BZ.
        """

        NKPTS = len(kpoints)

        # self.SW = np.array([self.spectral_weight_k(kpoints[ik])
        #                     for ik in range(NKPTS)], dtype=float)
        sw = []
        for ispin in range(self.wfc._nspin):
            if self.wfc._nspin == 2:
                if self.verbose:
                    print('#' * 60)
                    print('Spin component: %d' % ispin)
                    print('#' * 60)
            sw.append([self.spectral_weight_k(kpoints[ik], whichspin=ispin + 1) for ik in range(NKPTS)])

        self.SW = np.array(sw)

        # For noncollinear calculation, nspin = 1.
        # if self._lsoc:
        #     # self.SW = np.swapaxes(self.SW, 0, 1)
        #     self.SW = np.array([self.SW[0,:,:,0,:], self.SW[0,:,:,1,]])

        return self.SW

    def spectral_function(self, nedos=4000, sigma=0.02):
        r"""
        Generate the spectral function

            A(k_i, E) = \sum_m P_{Km}(k_i)\Delta(E - Em)

        Where the \Delta function can be approximated by Lorentzian or Gaussian
        function.
        """

        assert self.SW is not None, 'Spectral weight must be calculated first!'

        NS = self.wfc._nspin
        # Number of kpoints
        nk = self.SW.shape[1]
        # spectral function
        SF = np.zeros((NS, nedos, nk), dtype=float)

        emin = self.SW[:, :, :, 0].min()
        emax = self.SW[:, :, :, 0].max()
        e0 = np.linspace(emin - 5 * sigma, emax + 5 * sigma, nedos)

        for ispin in range(NS):
            for ii in range(nk):
                E_Km = self.SW[ispin, ii, :, 0]
                P_Km = self.SW[ispin, ii, :, 1]

                SF[ispin, :, ii] = np.sum(LorentzSmearing(e0[:, np.newaxis], E_Km[np.newaxis, :], sigma=sigma) * P_Km[np.newaxis, :],
                                          axis=1)
        return e0, SF


def spectral_weight_multiple_source(kpoints, unfold_objs, transform_matrix):
    """
    Calculate the spectral weight for a list of kpoints in the primitive BZ
    from a list of WAVECARs.
    """

    nk = len(kpoints)
    ns = unfold_objs[0].wfc._nspin
    for obj in unfold_objs[1:]:
        assert ns == obj.wfc._nspin

    # When reading from multiple WAVECAR, it is possible that each of them may have differnt number of bands
    # if so, we take only the first N bands, where N is the minimum values of bands
    nb = []
    for source in unfold_objs:
        nb.append(source.bands.shape[2])
    nbands = min(nb)

    spectral_weights = []
    for ispin in range(ns):
        sw_this_spin = []
        for ik in range(nk):
            this_k = kpoints[ik]
            this_k_supercell, _ = find_K_from_k(this_k, transform_matrix)
            # find index of K0
            ok = False
            for source in unfold_objs:
                try:
                    source.find_K_index(this_k_supercell)
                except ValueError:
                    continue
                sw_this_spin.append(source.spectral_weight_k(this_k, whichspin=ispin + 1)[:nbands, :])
                ok = True
                break

            if not ok:
                raise ValueError(f'This kpoint {this_k} this not found in the WAVECAR(s) provided')

        spectral_weights.append(sw_this_spin)

    return np.array(spectral_weights)


def concatenate_scf_kpoints(scf_kpts, scf_weights, kpoints):
    """Concatenate SCF kpoints (from IBZKPT) with zero-weighted kpoints"""
    scf_kpts = np.asarray(scf_kpts)
    kpoints = np.concatenate([scf_kpts, kpoints], axis=0)
    weights = np.zeros(kpoints.shape[0])
    weights[:len(scf_weights)] = scf_weights
    return kpoints, weights
