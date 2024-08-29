"""
Plotting utilities
"""
import os
from typing import Union, Sequence
import warnings
import itertools
import colorsys
import matplotlib.pyplot as plt

import numpy as np
from colormath.color_conversions import convert_color
from colormath.color_objects import (
    HSVColor,
    LabColor,
    LCHabColor,
    LCHuvColor,
    XYZColor,
    sRGBColor,
)
from matplotlib import cycler, rcParams
from matplotlib.style import context
from matplotlib.colors import to_rgb, cnames
from matplotlib.patches import Patch

from .unfold import UnfoldKSet, clean_latex_string, process_projection_options, parse_atoms
from .effective_mass import EffectiveMass, fitted_band

# pylint: disable=too-many-locals, too-many-arguments, import-outside-toplevel, too-many-branches, too-many-statements, too-many-nested-blocks, unexpected-keyword-arg


class UnfoldPlotter:
    """A collection plotting tools for unfolded band structures"""

    def __init__(self, unfold: UnfoldKSet):
        """
        Instantiate the plotter object

        :param unfold: `UnfoldKSet` object to be *plotted*.
        """
        self.unfold = unfold

    @staticmethod
    def plot_dos(ax, dos_plotter, dos_label, dos_options, ylim, eref, atoms=None, colours=None, orbitals_subplots=None):
        """
        Prepare and plot the density of states.
        """
        from pymatgen.electronic_structure.core import Spin

        if not dos_options:
            dos_options = {}

        from sumo.plotting import sumo_base_style, sumo_bs_style
        plt.style.use([sumo_base_style, sumo_bs_style])

        if atoms:
            # set DOS element & orbital colours to match the easyunfold band structure projections
            def _get_orbital_colour_dict(index, colour_list):
                # set s,p,d,f to different shades of colours[i]
                return {
                    's': adjust_lightness(colour_list[index], 1.0),
                    'p': adjust_lightness(colour_list[index], 0.7),
                    'px': adjust_lightness(colour_list[index], 0.7),
                    'py': adjust_lightness(colour_list[index], 0.8),
                    'pz': adjust_lightness(colour_list[index], 0.6),
                    'd': adjust_lightness(colour_list[index], 0.45),
                    'dxy': adjust_lightness(colour_list[index], 0.45),
                    'dyz': adjust_lightness(colour_list[index], 0.55),
                    'dxz': adjust_lightness(colour_list[index], 0.35),
                    'dz2': adjust_lightness(colour_list[index], 0.65),
                    'dx2-y2': adjust_lightness(colour_list[index], 0.25),
                    'x2-y2': adjust_lightness(colour_list[index], 0.25),  # labelled differently in VASP PROCAR
                    'f': adjust_lightness(colour_list[index], 0.2)
                }

            # Create a dictionary with different shades of colours[i] for each orbital
            sumo_colours = {atom: _get_orbital_colour_dict(i, colours) for i, atom in enumerate(atoms)}
            if len(atoms) != len(set(atoms)):
                # if atom entries are not all unique (i.e. repeated with different orbitals), then adjust colours of specific orbitals
                for i, atom, this_orbitals in zip(range(len(atoms)), atoms, orbitals_subplots):
                    if this_orbitals != 'all':
                        orbital_list = [token.strip() for token in this_orbitals.split(',')]
                        orbital_colour_dict = _get_orbital_colour_dict(i, colours)
                        for orbital in orbital_list:
                            # set all orbital keys with "orbital" in key:
                            for key in sumo_colours[atom].keys():
                                if orbital in key:
                                    sumo_colours[atom][key] = orbital_colour_dict[key]
        else:
            from pkg_resources import Requirement, resource_filename
            try:
                import configparser
            except ImportError:
                import ConfigParser as configparser

            config_path = resource_filename(Requirement.parse('sumo'), 'sumo/plotting/orbital_colours.conf')
            sumo_colours = configparser.ConfigParser()
            sumo_colours.read(os.path.abspath(config_path))

        # don't use first 4 colours; these are the band structure line colours:
        cycle = cycler('color', rcParams['axes.prop_cycle'].by_key()['color'][4:])
        with context({'axes.prop_cycle': cycle}):
            try:
                plot_data = dos_plotter.dos_plot_data(xmin=ylim[0],
                                                      xmax=ylim[1],
                                                      zero_energy=eref,
                                                      zero_to_efermi=False,
                                                      colours=sumo_colours,
                                                      **dos_options)
            except TypeError:  # sumo < 2.3
                try:
                    plot_data = dos_plotter.dos_plot_data(xmin=ylim[0],
                                                          xmax=ylim[1],
                                                          ref_energy=eref,
                                                          zero_to_efermi=False,
                                                          colours=sumo_colours,
                                                          **dos_options)
                except TypeError:  # sumo < 2.2
                    plot_data = dos_plotter.dos_plot_data(xmin=ylim[0],
                                                          xmax=ylim[1],
                                                          zero_to_efermi=True,
                                                          colours=sumo_colours,
                                                          **dos_options)

        mask = plot_data['mask']
        energies = plot_data['energies'][mask]
        lines = plot_data['lines']
        spins = [Spin.up] if len(lines[0][0]['dens']) == 1 else [Spin.up, Spin.down]

        for line_set in plot_data['lines']:
            for line, spin in itertools.product(line_set, spins):
                if spin == Spin.up or len(spins) == 1:
                    label = line['label']
                    densities = line['dens'][spin][mask]
                else:
                    label = ''
                    densities = -line['dens'][spin][mask]

                line_style = '-'
                if label:
                    if any(i in label for i in ['py', 'dyz']):
                        line_style = '--'
                    elif any(i in label for i in ['pz', 'dxz']):
                        line_style = ':'
                    elif any(i in label for i in ['dz2']):
                        line_style = '-.'
                    elif any(i in label for i in ['x2-y2']):
                        line_style = (5, (10, 3))

                ax.fill_betweenx(energies, densities, 0, lw=0, facecolor=line['colour'], alpha=line['alpha'])
                ax.plot(densities, energies, label=label, color=line['colour'], ls=line_style, lw=1)

        # x and y axis reversed versus normal dos plotting
        ax.set_ylim(*ylim)
        if len(spins) == 1:
            ax.set_xlim(0, plot_data['ymax'])
        else:
            ax.set_xlim(plot_data['ymin'], plot_data['ymax'])

        if dos_label is not None:
            ax.set_xlabel(dos_label)

        ax.set_yticks([])  # no y ticks
        ax.set_xticks([])  # no x ticks
        ax.legend(loc=2, frameon=False, ncol=1, bbox_to_anchor=(1.0, 1.0), fontsize=9)

        return ax

    def plot_spectral_function(
        self,
        engs: np.ndarray,
        sf: np.ndarray,
        dos_plotter=None,
        dos_label=None,
        dos_options=None,
        zero_line=False,
        eref=None,
        figsize=(4, 3),
        ylim=(-5, 5),
        dpi=300,
        vscale=1.0,
        contour_plot=False,
        alpha=1.0,
        save=False,
        ax=None,
        vmin=None,
        vmax=None,
        cmap='PuRd',
        show=False,
        title=None,
        intensity=1.0,
    ):
        """
        Plot the spectral function.

        :param engs: The energies of the spectral functions.
        :param sf: An array of the spectral function.
        :param eref: Reference energy to be used - this energy will be set as zero.
        :param figsize: Size of the figure.
        :param ylim: Plotting limit for the y-axis, with respect to `eref`.
        :param dpi: DPI of the generated graph.
        :param vscale: A normalisation/scaling factor for the colour map. Smaller values will increase the colour intensity.
        :param contour_plot: Whether to use contour plot instead of normal meshed color map.
        :param alpha: Color map transparency factor.
        :param save: The file to save the generated figure to.
        :param ax: An existing plotting axis to be used.
        :param vmin: Min value for the colour map.
        :param vmax: Max value for the colour map.
        :param cmap: Name of the color map to be used.

        :returns: The figure generated containing the spectral function.
        """
        unfold = self.unfold
        kdist = unfold.get_kpoint_distances()
        # Extend of x axis
        xmin, xmax = kdist.min(), kdist.max()
        nspin = sf.shape[0]

        if eref is None:
            eref = unfold.calculated_quantities.get('vbm', 0.0)

        if dos_plotter:
            fig, axes = plt.subplots(1, 2, facecolor='w', gridspec_kw={'width_ratios': [3, 1], 'wspace': 0}, figsize=figsize, dpi=dpi)
            if nspin > 1:
                warnings.warn('DOS plotter is not supported for spin-separated plots. Reverting to non spin-polarised plotting.')
                nspin = 1
        elif ax is None:
            if nspin == 1:
                fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
                axes = [ax]
            else:
                fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
        else:
            axes = [ax] if not isinstance(ax, list) else ax
            fig = axes[0].figure

        # Shift the kdist so the pcolormesh draw the pixel centred on the original point
        X, Y = np.meshgrid(kdist, engs - eref)

        # Calculate the min and max values within the field of view, scaled by the factor
        mask = (engs < (ylim[1] + eref)) & (engs > (ylim[0] + eref))
        vmin = sf[:, mask, :].min() if vmin is None else vmin

        if vmax is None:
            vmax = sf[:, mask, :].max()
            vmax = (vmax - vmin) * (vscale /
                                    intensity) + vmin  # vscale and intensity have the equal opposite effect on the colour intensity

        for ispin, ax_ in zip(range(nspin), axes):
            if contour_plot:
                ax_.contourf(X, Y, sf[ispin], cmap=cmap, vmax=vmax, vmin=vmin, alpha=alpha)
            else:
                ax_.pcolormesh(X, Y, sf[ispin], cmap=cmap, shading='auto', vmax=vmax, vmin=vmin, alpha=alpha)

            ax_.set_xlim(xmin, xmax)
            ax_.set_ylim(*ylim)
            ax_.set_ylabel('Energy (eV)', labelpad=5)
            if title:
                ax_.set_title(title)

            # Label the kpoints
            self._add_kpoint_labels(ax_)

        if dos_plotter:
            ax = fig.axes[1]
            ax = self.plot_dos(ax, dos_plotter, dos_label, dos_options, ylim, eref)

        if zero_line:
            try:
                from sumo.plotting import draw_themed_line
                try:
                    draw_themed_line(0, axes[0], zorder=10)  # bump zorder to put on top of pcolormesh
                except TypeError:  # old sumo
                    draw_themed_line(0, axes[0])
                if dos_plotter:
                    draw_themed_line(0, axes[1])

            except ImportError:
                warnings.warn('zero_line option requires sumo to be installed!')

        fig.tight_layout(pad=0.2)
        if save:
            fig.savefig(save, dpi=dpi, bbox_inches='tight')
        if show:
            fig.show()
        return fig

    def _plot_spectral_function_rgba(
        self,
        engs: np.ndarray,
        sf: np.ndarray,
        eref: Union[None, float] = None,
        figsize=(4, 3),
        ylim=(-3, 3),
        dpi: float = 150,
        vscale: float = 1.0,
        save: bool = False,
        ax: Union[None, plt.Axes] = None,
        show: bool = False,
        title: Union[None, str] = None,
        vmin: Union[None, float] = None,
        vmax: Union[None, float] = None,
        intensity: float = 1.0,
    ):
        """
        Plot spectral function defined as RGBA colours

        :param np.ndarray engs: The energies of the spectral functions.
        :param sf: An array of the spectral function.
        :param float eref: Reference energy to be used - this energy will be set as zero.
        :param figsize: Size of the figure.
        :param ylim: Plotting limit for the y-axis, with respect to `eref`.
        :param dpi: DPI of the generated graph.
        :param save: The file where the generated figure is saved.
        :param ax: An existing plotting axis to be used.
        :param vmin: Lower bound for constructing the alpha channel.
        :param vmax: Upper bound for constructing the alpha channel.
        :param vscale: Scaling factor for the alpha channel.
        :param title: Title for the plot.

        :returns: the figure generated containing the spectral function.
        """
        unfold = self.unfold
        nspin = sf.shape[0]

        if eref is None:
            eref = unfold.calculated_quantities.get('vbm', 0.0)

        if ax is None:
            if nspin == 1:
                fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
                axes = [ax]
            else:
                fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
        else:
            axes = ax if isinstance(ax, (list, np.ndarray)) else [ax]
            fig = axes[0].figure

        mask = (engs < (ylim[1] + eref)) & (engs > (ylim[0] + eref))
        vmin = sf[:, mask, :, 3].min() if vmin is None else vmin

        if vmax is None:
            vmax = sf[:, mask, :, 3].max()

        # Clip the alpha range
        alpha = sf[:, :, :, 3]
        alpha = ((alpha - vmin) / (vmax - vmin)) * (intensity / vscale)  # vscale & intensity have equal opposite effects
        alpha = np.clip(alpha, 0, 1)
        sf[:, :, :, 3] = alpha

        for ispin, ax_ in zip(range(nspin), axes):
            # plt.imshow's pixel coordinates are not at the centre of the pixel
            # Hence, the `extent`` needs to be updated so the centre of the pixel is aligned with the coordinates:
            extent = np.array([0, sf.shape[2], max(engs) - eref, min(engs) - eref])
            ebin = (max(engs) - min(engs)) / sf.shape[1]
            extent[:2] -= 0.5
            extent[2:] -= ebin * 0.5
            ax_.imshow(sf[ispin], extent=extent, aspect='auto', origin='upper')
            ax_.set_ylim(ylim)
            ax_.set_xlim(0, sf.shape[2])
            ax_.set_ylabel('Energy (eV)', labelpad=5)
            ax_.set_title(title)
            self._add_kpoint_labels(ax_, x_is_kidx=True)

        fig.tight_layout(pad=0.2)
        if save:
            fig.savefig(save, dpi=dpi)
        if show:
            fig.show()
        return fig

    def _add_kpoint_labels(self, ax: plt.Axes, x_is_kidx=False):
        """Add labels to the k-points for a given axes"""
        # Label the kpoints
        labels = self.unfold.get_combined_kpoint_labels()
        kdist = self.unfold.get_kpoint_distances()

        # Explicit label indices
        tick_locs = []
        tick_labels = []
        for index, label in labels:
            xloc = index if x_is_kidx else kdist[index]
            ax.axvline(x=xloc, lw=0.5, color='k', ls=':', alpha=0.8)
            tick_locs.append(xloc)
            tick_labels.append(clean_latex_string(label))
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labels)

    def plot_effective_mass(self,
                            eff: EffectiveMass,
                            engs: np.ndarray,
                            sf: np.ndarray,
                            eref: Union[None, float] = None,
                            save: Union[None, str] = None,
                            show: bool = False,
                            effective_mass_data: dict = None,
                            **kwargs):
        """
        Plot the effective masses on top of the spectral function.

        :param eff: An `EffectiveMass` object used for plotting.
        :param engs: The energies of the spectral functions.
        :param sf: An array of the spectral function.
        :param eref: Reference energy to be used - this energy will be set as zero.
        :param effective_mass_data: Calculated data of the effective masses.

        :returns: A figure with the data used for fitting effective mass plotted on top of the spectral function.
        """

        kcbm = eff.get_band_extrema(mode='cbm')[0]
        kvbm = eff.get_band_extrema(mode='vbm')[0]
        all_k = sorted(list(set(list(kcbm) + list(kvbm))))
        kdist = self.unfold.get_kpoint_distances()
        if eref is None:
            eref = self.unfold.calculated_quantities['vbm']

        fig, axes = plt.subplots(1, len(all_k), figsize=(4 * len(all_k), 3), dpi=300, squeeze=False)
        if effective_mass_data is None:
            effective_mass_data = eff.get_effective_masses()

        # Plot the spectral function
        xwidth = abs(kdist[1] - kdist[0])
        for (ik, ax) in zip(all_k, axes[0]):
            self.plot_spectral_function(engs, sf, ax=ax, eref=eref, **kwargs)
            xk = kdist[ik]
            xlim = (xk - xwidth / 2, xk + xwidth / 2)
            ax.set_xlim(xlim)
            ax.set_title(f'Kpoint: {ik}')

        elec = effective_mass_data.get('electrons', [])
        # Plot the detected effective mass fitting data on top
        for entry in elec:
            ik = entry['kpoint_index']
            iax = all_k.index(ik)
            x = entry['raw_data']['raw_fit_values'][0]
            y = entry['raw_data']['raw_fit_values'][1]
            axes[0, iax].plot(x, np.asarray(y) - eref, '-o', color='C1')
            axes[0, iax].set_xlim(min(x) - xwidth / 2, max(x) + xwidth / 2)

        hole = effective_mass_data.get('holes', [])
        for entry in hole:
            ik = entry['kpoint_index']
            iax = all_k.index(ik)
            x = entry['raw_data']['raw_fit_values'][0]
            y = entry['raw_data']['raw_fit_values'][1]
            axes[0, iax].plot(x, np.asarray(y) - eref, '-o', color='C2')
            axes[0, iax].set_xlim(min(x) - xwidth / 2, max(x) + xwidth / 2)
        if save:
            fig.savefig(save)

        if show:
            fig.show()
        return fig

    def plot_spectral_weights(
        self,
        figsize=(4, 3),
        ylim=(-3, 3),
        dpi: float = 150,
        factor: float = 3.0,
        eref: Union[None, float] = None,
        color: str = 'C1',
        alpha: float = 0.5,
        save: Union[None, str] = None,
        ax: Union[None, plt.Axes] = None,
        show: bool = False,
        title: Union[None, str] = None,
    ):
        """
        Plot the spectral weights.

        ```{note}
        The reduction of symmetry means there can be multiple supercell kpoints for each
        primitive cell kpoint. When using this scattering plot representation, the markers can
        overlap each other leading to misrepresentations of the actual effective band structure.
        ```

        However, this function is still useful when: 1. the symmetry splitting is turned off. 2.
        direct visualization of the underlying spectral weight is needed. 3. Check the correctness
        of effective mass extraction.

        :param eref: Reference energy to be used - this energy will be set as zero.
        :param figsize: Size of the figure.
        :param ylim: Plotting limit for the y-axis, with respect to `eref`.
        :param dpi: DPI of the generated graph.
        :param alpha: Alpha for the markers.
        :param save: Name of the file where the generated figure is saved.
        :param ax: Existing plotting axes to be used (list if having two spin channels).
        :param factor: Scaling factor for the marker size.
        :param color: Color of the markers.
        :param title: Title for the generated plot.

        :returns: A Figure with the spectral weights plotted as a scatter plot.
        """
        unfold = self.unfold
        kdist = unfold.get_kpoint_distances()
        # a list of spectral weights each with (nspin, nkset, nbands, 2)
        sws = unfold.get_spectral_weights()
        nspin = sws[0].shape[0]

        nb = sws[0].shape[2]

        if eref is None:
            eref = unfold.calculated_quantities.get('vbm', 0.0)

        if ax is None:
            if nspin == 1:
                fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
                axes = [ax]
            else:
                fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
        else:
            axes = [ax] if not isinstance(ax, list) else ax
            fig = axes[0].figure

        kweights = unfold.expansion_results['weights']
        xs = []
        ys = []
        ss = []
        # Each spin channel
        for (ispin, ax_) in zip(range(nspin), axes):
            # Collect data
            for (ik, sw) in enumerate(sws):
                # Each sub kpoint
                for isubk in range(sw.shape[1]):
                    x = np.repeat(kdist[ik], nb)
                    y = sw[ispin, isubk, :, 0] - eref
                    marker_size = (sw[ispin, isubk, :, 1] * factor) * kweights[ik][isubk]
                    xs.append(x)
                    ys.append(y)
                    ss.append(marker_size)

            ax_.scatter(
                np.vstack(xs),
                np.vstack(ys),
                s=np.vstack(ss),
                color=color,
                alpha=alpha,
                lw=0.0,
            )
            ax_.set_xlim(0, kdist.max())
            ax_.set_ylim(ylim)
            ax_.set_ylabel('Energy [eV]', labelpad=5)
            if title:
                ax_.set_title(title)
            self._add_kpoint_labels(ax_)

        fig.tight_layout(pad=0.2)

        if save:
            fig.savefig(save, dpi=dpi)
        if show:
            fig.show()
        return fig

    def plot_projected(
        self,
        procar: Union[str, list] = 'PROCAR',
        dos_plotter=None,
        dos_label=None,
        dos_options=None,
        zero_line=False,
        eref=None,
        gamma=False,
        npoints=2000,
        sigma=0.2,
        ncl=False,
        symm_average=True,
        figsize=(4, 3),
        ylim=(-5, 5),
        dpi=300,
        vscale=1.0,
        contour_plot=False,
        alpha=1.0,
        save=False,
        ax=None,
        cmap='PuRd',
        show=False,
        title=None,
        atoms=None,
        poscar='POSCAR',
        atoms_idx=None,
        orbitals=None,
        use_subplot=False,
        colours=None,
        colorspace='lab',
        intensity=1.0,
    ):
        """
        Plot projected spectral function onto multiple subplots or a single plot with color mapping.

        This simply computes the spectral function at each orbital/atoms sites and plot them onto
        multiple subplots. The columns are for each orbital and the rows are for each spin channel.

        :param procar: Name(s) of the `PROCAR(.gz)` file(s).

        :param colours: Default is pastel red, green, blue if <=3 projections, else red, green,
            blue, purple, orange, yellow.

        :returns: Generated plot.
        """
        unfoldset = self.unfold
        unfoldset.load_procars(procar)
        nspin = unfoldset.calculated_quantities['spectral_weights_per_set'][0].shape[0]

        if atoms_idx is not None:
            atoms_idx_subplots = atoms_idx.split('|')  # list of strings

        else:  # parse atoms
            atoms, atoms_idx_subplots, _ = parse_atoms(atoms, orbitals, poscar)
            # atoms as a list, atoms_idx_subplots as a list of lists

        if orbitals is not None:
            orbitals_subplots = orbitals.split('|')

            # Special case: if only one set is passed, apply it to all atomic specifications
            if len(orbitals_subplots) == 1:
                orbitals_subplots = orbitals_subplots * len(atoms_idx_subplots)

            if len(orbitals_subplots) != len(atoms_idx_subplots):
                raise ValueError('The number of elected orbitals and atoms indices are not matched.')
        # If not set, use all for all subsets
        else:
            orbitals_subplots = ['all'] * len(atoms_idx_subplots)

        # Load the data

        nsub = len(atoms_idx_subplots)

        all_sf = []
        vmaxs = []

        if eref is None:
            eref = unfoldset.calculated_quantities.get('vbm', 0.0)

        # Collect spectral functions and scale
        for this_idx, this_orbitals in zip(atoms_idx_subplots, orbitals_subplots):
            # Set up the atoms_idx and orbitals
            if isinstance(this_idx, str):
                this_idx, this_orbitals = process_projection_options(this_idx, this_orbitals)
            elif this_orbitals != 'all':
                this_orbitals = [token.strip() for token in this_orbitals.split(',')]

            eng, spectral_function = unfoldset.get_spectral_function(gamma=gamma,
                                                                     npoints=npoints,
                                                                     sigma=sigma,
                                                                     ncl=ncl,
                                                                     atoms_idx=this_idx,
                                                                     orbitals=this_orbitals,
                                                                     symm_average=symm_average)
            all_sf.append(spectral_function)

            emin, emax = ylim
            # Clip the effective range
            mask = (eng < (emax + eref)) & (eng > (emin + eref))
            vmin = spectral_function[:, mask, :].min()
            vmax = spectral_function[:, mask, :].max()
            vmax = (vmax - vmin) * (vscale / intensity) + vmin  # vscale & intensity have equal opposite effects
            vmaxs.append(vmax)

        # Workout the vmax and vmin
        vmax = max(vmaxs)
        vmin = 0.

        if use_subplot:
            if dos_plotter:
                warnings.warn('DOS plotter is not supported for projected subplots. Use `--combined` if '
                              'you want to plot the DOS along with the projected band structure!')
            fig, axs = plt.subplots(nspin, len(atoms_idx_subplots), sharex=True, sharey=True, squeeze=False, figsize=(3.0 * nsub, 4.0))
            # Plot the spectral function with constant colour scales
            for spectral_function, i in zip(all_sf, range(len(atoms_idx_subplots))):
                if (title is None or any(atom == title for atom in atoms)) and atoms is not None:
                    # use atoms as titles (checks if atom == title from previous subplot)
                    title = atoms[i]

                plotter = UnfoldPlotter(unfoldset)
                plotter.plot_spectral_function(
                    eng,
                    spectral_function,
                    eref=eref,
                    save=save,
                    vscale=vscale,
                    vmax=vmax,
                    vmin=vmin,
                    cmap=cmap,
                    title=title,
                    ax=axs[:, i].tolist(),
                    figsize=figsize,
                    show=show,
                    contour_plot=contour_plot,
                    alpha=alpha,
                    dpi=dpi,
                    ylim=ylim,
                    intensity=intensity,
                )
        else:
            # Make a combined plot
            sf_size = all_sf[0].shape
            stacked_sf = np.stack(all_sf, axis=-1).reshape(np.prod(sf_size), len(all_sf))

            # Construct the colour basis
            if colours is None:
                if len(all_sf) <= 3:
                    colours = ['#CC33A7', '#A7CC33', '#33A7CC']
                else:
                    colours = [
                        (1, 0, 0),  # red
                        (0, 1, 0),  # green
                        (0, 0, 1),  # blue
                        (152 / 255, 78 / 255, 163 / 255),  # purple
                        (1, 127 / 255, 0),  # orange
                        (1, 1, 51 / 255),  # yellow
                    ]
            colours = colours[:len(all_sf)]
            # Compute spectral weight data with RGB reshape it back into the shape (nengs, nk, 3)
            sf_rgb = interpolate_colors(colours, stacked_sf, colorspace, normalize=True).reshape(sf_size + (3,))
            sf_sum = np.sum(all_sf, axis=0)[:, :, :, None]
            sf_rgba = np.concatenate([sf_rgb, sf_sum], axis=-1)

            if dos_plotter:
                fig, axes = plt.subplots(1,
                                         2,
                                         facecolor='w',
                                         gridspec_kw={
                                             'width_ratios': [3, 1],
                                             'wspace': 0
                                         },
                                         figsize=figsize,
                                         dpi=dpi,
                                         squeeze=True)
                ax = axes[0]
                if nspin > 1:
                    warnings.warn('DOS plotter is not supported for spin-separated plots. Reverting to non spin-polarised plotting.')
                    nspin = 1

            elif ax is None:
                fig, ax = plt.subplots(1, nspin, figsize=figsize, squeeze=True)

            self._plot_spectral_function_rgba(
                eng,
                sf_rgba,
                eref=eref,
                save=save,
                title=title,
                vscale=vscale,
                ax=ax,
                figsize=figsize,
                show=show,
                dpi=dpi,
                ylim=ylim,
                intensity=intensity,
            )

            if dos_plotter:
                ax = fig.axes[1]
                ax = self.plot_dos(ax, dos_plotter, dos_label, dos_options, ylim, eref, atoms, colours, orbitals_subplots)

            if zero_line:
                try:
                    from sumo.plotting import draw_themed_line
                    if dos_plotter:
                        try:
                            draw_themed_line(0, axes[0], zorder=10)  # bump zorder to put on top of pcolormesh
                        except TypeError:  # old sumo
                            draw_themed_line(0, axes[0])
                        draw_themed_line(0, axes[1])
                    else:
                        draw_themed_line(0, ax)

                except ImportError:
                    warnings.warn('zero_line option requires sumo to be installed!')

            if atoms is not None:  # add figure legend with atoms and colors
                legend_elements = [Patch(facecolor=colours[i], label=atom, alpha=0.7) for i, atom in enumerate(atoms)]
                fig.axes[0].legend(handles=legend_elements, bbox_to_anchor=(1.025, 1), fontsize=9)
                fig.subplots_adjust(right=0.78)  # ensure legend is not cut off

        return fig

    @staticmethod
    def plot_effective_mass_fit(efm: EffectiveMass,
                                npoints: int = 3,
                                carrier: str = 'electrons',
                                idx: int = 0,
                                ax: Union[plt.Axes, None] = None,
                                save: Union[None, str] = None,
                                dpi: float = 120):
        """
        Plot detected band edges and the fitted effective masses.

        :param efm: `EffectiveMass` object for plotting.
        :param npoints: The number of points to be used for fitting.
        :param carrier: Type of the charge carrier, e.g. `electrons` or `holes`.
        :param idx: Index for the detected effective mass of the same kind.
        :param ax: `plt.Axes` used for plotting.
        :param save: Name of the file used for saveing.
        :param dpi: DPI of the figure when saving.

        :returns: A figure with plotted data.
        """
        data = efm.get_effective_masses(npoints=npoints)[carrier][idx]

        x = data['raw_data']['kpoint_distances']
        y = data['raw_data']['effective_energies']
        me = data['effective_mass']
        y1 = fitted_band(x, me)
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.plot(x, y, 'x-', label='Energy ')
        ax.plot(x, y1, label='fitted')
        ax.legend()
        if save:
            fig.savefig(save, dpi=dpi)
        return fig


def interpolate_colors(colours: Sequence, weights: list, colorspace='lab', normalize=True):
    """
    Interpolate colours at a number of points within a colorspace.

    :param colours: A list of colours specified in any way supported by matplotlib.
    :param weights: A list of weights with the shape (n, N). Where the N values of
      the last axis give the amount of N colours supplied in `colours`.
    :param colorspace: The colorspace in which to perform the interpolation. The
            allowed values are rgb, hsv, lab, luvlch, lablch, and xyz.

    :returns: A list of colours, specified in the rgb format as a (n, 3) array.
    """

    # Set up and check the color spaces
    colorspace_mapping = {
        'rgb': sRGBColor,
        'hsv': HSVColor,
        'lab': LabColor,
        'luvlch': LCHuvColor,
        'lablch': LCHabColor,
        'xyz': XYZColor,
    }

    if colorspace not in colorspace_mapping:
        raise ValueError(f'colorspace must be one of {colorspace_mapping.keys()}')

    colorspace = colorspace_mapping[colorspace]

    # Convert matplotlib color specification to colormath sRGB
    colors_srgb = [sRGBColor(*to_rgb(c)) for c in colours]

    colors_basis = [np.array(convert_color(srgb, colorspace, target_illuminant='d50').get_value_tuple()) for srgb in colors_srgb]

    # ensure weights is a numpy array
    weights = np.asarray(weights)

    # Normalise the weights if needed
    if normalize:
        weights = weights / np.sum(weights, axis=1)[:, None]  # each row sums to 1

    # perform the interpolation in the colorspace basis
    interpolated_colors = colors_basis[0] * weights[:, 0][:, None]
    for i in range(1, len(colors_basis)):
        interpolated_colors += colors_basis[i] * weights[:, i][:, None]

    # convert the interpolated colors back to RGB
    rgb_colors = [convert_color(colorspace(*c), sRGBColor).get_value_tuple() for c in interpolated_colors]

    # ensure all rgb values are less than 1 (sometimes issues in interpolation)
    normalised_rgb_colors = []
    for rgb_color_tuple in rgb_colors:
        if np.max(rgb_color_tuple) > 1:
            normalised_rgb_color = np.array(rgb_color_tuple) / np.max(rgb_color_tuple)
        else:
            normalised_rgb_color = np.array(rgb_color_tuple)

        normalised_rgb_color = np.clip(normalised_rgb_color, 0, 1)  # ensure all rgb values are between 0 and 1
        # if too white, darken:
        if np.linalg.norm(normalised_rgb_color) > 1:  # white af
            normalised_rgb_color *= (1 / np.linalg.norm(normalised_rgb_color)**(1 / 2))

        normalised_rgb_colors.append(normalised_rgb_color)

    rgb_colors = np.stack(normalised_rgb_colors, axis=0)

    return rgb_colors


def adjust_lightness(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> adjust_lightness('g', 1.1)  # slightly lighter
    >> adjust_lightness('r', 0.6)  # significantly darker
    >> adjust_lightness('#F034A3', 0.6)  # slightly darker
    >> adjust_lightness((.3,.55,.1), 1.2)  # slightly lighter
    """

    try:
        col = cnames[color]
    except KeyError:
        col = color
    col = colorsys.rgb_to_hls(*to_rgb(col))
    return colorsys.hls_to_rgb(col[0], 1 - amount * (1 - col[1]), col[2])
