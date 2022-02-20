#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The main module for unfolding workflow and algorithm
"""
# pylint: disable=invalid-name,protected-access

############################################################
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import spglib
from monty.json import MSONable

from easyunfold import __version__
from .pyvaspwfc.vaspwfc import vaspwfc

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
    out_weights = np.array(out_weights) / sum(out_weights)
    return out_points, out_weights


class UnfoldKSet(MSONable):
    """Stores the information of the kpoints in the primitive cell, and what they unfolds to in the supercell"""

    def __init__(self,
                 M,
                 kpts_pc,
                 pc_latt,
                 pc_opts,
                 sc_opts,
                 time_reversal=True,
                 expansion_results=None,
                 calculated_quantities=None,
                 kpoint_labels=None):
        """
        Args:
            kpts_pc: A list of kpoints in the PC
            pc_lattice: A 3x3 matrix of row lattice vectors of the primitive cell
            pc_opts: Symmetry operations of the primitive cell
            sc_opts: Symmetry operations of the supercell
        """
        # Basic properties - needed to recreate the object
        self.kpts_pc = kpts_pc
        self.pc_latt = pc_latt
        self.pc_opts = pc_opts
        self.sc_opts = sc_opts
        self.M = M
        self.expansion_results = expansion_results
        self.time_reversal = time_reversal
        self.calculated_quantities = {} if not calculated_quantities else calculated_quantities
        self.kpoint_labels = kpoint_labels
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
    def from_atoms(cls, M, kpts_pc, pc, sc, time_reversal=True):
        """Initialise from primitive cell and supercell atoms"""
        pc_symm_data = get_symmetry_dataset(pc)
        sc_symm_data = get_symmetry_dataset(sc)
        return cls(
            M=M,
            kpts_pc=kpts_pc,
            pc_latt=np.asarray(pc.cell),
            pc_opts=pc_symm_data['rotations'],
            sc_opts=sc_symm_data['rotations'],
            time_reversal=time_reversal,
        )

    def expand_pc_kpoints(self):
        """Comptue the pc kpoints to be unfolded into"""
        expended_k = []
        expended_weights = []
        for kpt in self.kpts_pc:
            kset, weights = expand_K_by_symmetry(kpt, self.pc_opts, self.sc_opts, time_reversal=self.time_reversal)
            expended_k.append(kset)
            expended_weights.append(weights)
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

        # Mapping between the pckpts to the redcued sckpts
        reduced_sc_map = []
        for sc_set in expended_sc:
            map_indx = []
            for _ in sc_set:
                map_indx.append(sc_kpts_map.pop(0))
            reduced_sc_map.append(map_indx)

        self.expansion_results['reduced_sckpts'] = reduced_sckpts
        self.expansion_results['reduced_sckpts_map'] = reduced_sc_map
        # A nested list that stores the indices of the sc kpts in the reduced_sckpts list

    def write_sc_kpoints(self, file):
        """Write the supercell kpoints"""
        if not self.expansion_results.get('reduced_sckpts'):
            self.generate_sc_kpoints()
        write_kpoints(self.expansion_results['reduced_sckpts'], file, comment='supercell kpoints')

    def write_pc_kpoints(self, file, expanded=False):
        """Write the primitive cell kpoints"""
        if expanded:
            all_pc = []
            for tmp in self.expansion_results['kpoints']:
                all_pc.extend(tmp)
            write_kpoints(all_pc, file, comment='expanded primitive cell kpoints')
        else:
            write_kpoints(self.kpts_pc, file, comment='primitive cell kpoints')

    def _get_spectral_weights(self,
                              wavecar,
                              npoints=2000,
                              sigma=0.1,
                              gamma=False,
                              lsorbit=False,
                              gamma_half='x',
                              also_spectral_function=False,
                              symm_average=True):
        """
        Fetch spectral weights from a wavecar and compute spectral function is requested
        """
        if wavecar is None:
            # Use existing results
            if self.calculated_quantities['spectral_weights_is_averaged'] != symm_average:
                tmp = 'not ' if not symm_average else ''
                tmp2 = bool(not symm_average)
                raise RuntimeError(f'Previously calculated spectral weights was {tmp}averaged. Please set symm_avg to be {tmp2}.')
            sws = self.calculated_quantities['spectral_weights']
        else:
            unfold_obj = Unfold(self.M, wavecar, gamma=gamma, lsorbit=lsorbit, gamma_half=gamma_half)
            sws = []
            if symm_average is True:
                for kset, weights in zip(self.expansion_results['kpoints'], self.expansion_results['weights']):
                    sw = unfold_obj.spectral_weight(kset)
                    for ik, w in enumerate(weights):
                        sw[:, ik, :, :] *= w
                    sws.append(sw.sum(axis=1))
                sws = np.stack(sws, axis=1)
            else:
                sws = unfold_obj.spectral_weight(self.kpts_pc)
            self.calculated_quantities['spectral_weights'] = sws
            self.calculated_quantities['spectral_weights_is_averaged'] = symm_average

        if also_spectral_function:
            e0, spectral_function = spectral_function_from_weights(sws, nedos=npoints, sigma=sigma)
            self.calculated_quantities['e0'] = e0
            self.calculated_quantities['spectral_function'] = spectral_function
            return sws, e0, spectral_function
        return sws

    def get_spectral_function(self, wavecar=None, npoints=2000, sigma=0.1, gamma=False, lsorbit=False, gamma_half='x', symm_average=True):
        """Get the spectral function"""
        _, e0, spectral_function = self._get_spectral_weights(wavecar,
                                                              npoints=npoints,
                                                              sigma=sigma,
                                                              gamma=gamma,
                                                              lsorbit=lsorbit,
                                                              gamma_half=gamma_half,
                                                              also_spectral_function=True,
                                                              symm_average=symm_average)
        return e0, spectral_function

    def get_spectral_weights(self, wavecar=None, gamma=False, lsorbit=False, gamma_half='x', symm_average=True):
        """Get the spectral function"""
        return self._get_spectral_weights(wavecar=wavecar,
                                          gamma=gamma,
                                          lsorbit=lsorbit,
                                          gamma_half=gamma_half,
                                          also_spectral_function=False,
                                          symm_average=symm_average)

    def as_dict(self):
        """To a dictionary representation"""
        output = {'@module': self.__class__.__module__, '@class': self.__class__.__name__, '@version': __version__}
        for key in [
                'M', 'kpts_pc', 'pc_latt', 'pc_opts', 'sc_opts', 'expansion_results', 'time_reversal', 'calculated_quantities',
                'kpoint_labels'
        ]:
            output[key] = getattr(self, key)
        return output


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


def write_kpoints(kpoints: np.ndarray, outpath='KPOINTS', comment=''):
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
    ik = 0
    for line in content[3:]:
        tokens = line.split()
        this_kpt = [float(value) for value in tokens[:3]]

        if len(tokens) >= 5:
            labels.append([ik, tokens[4]])

        kpts.append(this_kpt)
        ik += 1
        if ik == nkpts:
            break
    return kpts, comment, labels


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
                atomic_colors=None,
                eref=0.0,
                nseg=None,
                save='ebs_s.png',
                kpath_label=None,
                factor=20,
                figsize=(3.0, 4.0),
                ylim=(-3, 3),
                show=True,
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

    #mpl.rcParams['axes.unicode_minus'] = False
    atomic_colors = [] if atomic_colors is None else atomic_colors
    kpath_label = [] if kpath_label is None else kpath_label

    nspin = spectral_weight.shape[0]
    kpt_c = np.dot(kpts, np.linalg.inv(cell).T)
    kdist = np.r_[0, np.cumsum(np.linalg.norm(np.diff(kpt_c, axis=0), axis=1))]
    nb = spectral_weight.shape[2]
    # x0 = np.outer(np.ones(nb), kdist).T
    x0 = np.tile(kdist, (nb, 1)).T

    if atomic_weights is not None:
        atomic_weights = np.asarray(atomic_weights)
        assert atomic_weights.shape[1:] == spectral_weight.shape[:-1]

        if not atomic_colors:
            atomic_colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']

    fig = plt.figure()
    fig.set_size_inches(figsize)
    if nspin == 1:
        axes = [plt.subplot(111)]
        fig.set_size_inches(figsize)
    else:
        axes = [plt.subplot(121), plt.subplot(122)]
        fig.set_size_inches((figsize[0] * 2, figsize[1]))

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

    plt.tight_layout(pad=0.2)
    plt.savefig(save, dpi=360)
    if show:
        plt.show()


def EBS_cmaps(kpts,
              cell,
              E0,
              spectral_function,
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
              cmap='jet'):
    """
    plot the effective band structure with colormaps.  The plotting function
    utilizes Matplotlib package.

    Args:
        kpts: the kpoints vectors in fractional coordinates.
        cell: the primitive cell basis
        e0: The energies corresponds to each element of the spectral function
        spectral_function: The spectral function array in the shape of (nspin, nk, neng)
        eref: Refernce point for zero energy
        kpath_label: Label of the high symmetry kpoints along the pathway
        nseg: Number of points in each segment of the kpoint pathway
        explicit_labels: A list of tuplies containing tuples of `(index, label)` to explicitly set kpoint labels.
        save: Name of the file the plot to be saved to.
        figsize: SIze of hte figure
        ylim: Limit for the y axis. The limit is applied *after* substracting the refence energy.
        show: To show the plot interactively or not.
        contour_plot: Plot in the contour mode
        ax: Existing axis(axes) to plot onto
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
    for ispin in range(nspin):
        ax = axes[ispin]
        if contour_plot:
            ax.contourf(X, Y, spectral_function[ispin], cmap=cmap)
        else:
            ax.pcolormesh(X, Y, spectral_function[ispin], cmap=cmap, shading='auto')

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


def spectral_function_from_weights(spectral_weights, nedos=4000, sigma=0.02):
    r"""
    Generate the spectral function

        A(k_i, E) = \sum_m P_{Km}(k_i)\Delta(E - Em)

    Where the \Delta function can be approximated by Lorentzian or Gaussian
    function.
    """

    nk = spectral_weights.shape[1]
    ns = spectral_weights.shape[0]
    spectral_function = np.zeros((ns, nedos, nk), dtype=float)

    emin = spectral_weights[:, :, :, 0].min()
    emax = spectral_weights[:, :, :, 0].max()
    e0 = np.linspace(emin - 5 * sigma, emax + 5 * sigma, nedos)

    for ispin in range(ns):
        for ii in range(nk):
            E_Km = spectral_weights[ispin, ii, :, 0]
            P_Km = spectral_weights[ispin, ii, :, 1]

            spectral_function[ispin, :,
                              ii] = np.sum(LorentzSmearing(e0[:, np.newaxis], E_Km[np.newaxis, :], sigma=sigma) * P_Km[np.newaxis, :],
                                           axis=1)
    return e0, spectral_function


############################################################


class Unfold():
    """
    Unfold the band structure from Supercell calculation into a primitive cell and
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
        # all the KS energies
        self.bands = self.wfc._bands

        # G-vectors within the cutoff sphere, let's just do it once for all.
        # self.allGvecs = np.array([self.wfc.gvectors(ikpt=kpt+1)
        #                           for kpt in range(self.wfc._nkpts)], dtype=int)

        # spectral weight for all the kpoints
        self.SW = None

        self.verbose = verbose

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

    # def spectral_weight(self, kpoints, nproc=None):
    #     """
    #     Calculate the spectral weight for a list of kpoints in the primitive BZ.
    #     Here, we use "multiprocessing" package to parallel over the kpoints.
    #     """
    #
    #     NKPTS = len(kpoints)
    #
    #     if nproc is None:
    #         nproc = multiprocessing.cpu_count()
    #
    #     pool = multiprocessing.Pool(processes=nproc)
    #
    #     results = []
    #     for ik in range(NKPTS):
    #         res = pool.apply_async(self.spectral_weight_k, (kpoints[ik],))
    #         results.append(res)
    #
    #     self.SW = np.array([res.get() for res in results], dtype=float)
    #
    #     pool.close()
    #     pool.join()
    #
    #     return self.SW

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

        # NS = 2 if self._lsoc else self.wfc._nspin
        # For noncollinear calculation, nspin = 1.
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
