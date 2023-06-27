"""
Plotting utilities
"""
from typing import Union, Sequence
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
from matplotlib.colors import to_rgb

from .unfold import UnfoldKSet, clean_latex_string, process_projection_options
from .effective_mass import EffectiveMass, fitted_band

# pylint: disable=too-many-locals, too-many-arguments


class UnfoldPlotter:
    """A collection plotting tools for unfolded band structures"""

    def __init__(self, unfold: UnfoldKSet):
        """
        Instantiate the plotter object

        :param unfold: `UnfoldKSet` object to be *plotted*.
        """
        self.unfold = unfold

    def plot_spectral_function(
        self,
        engs: np.ndarray,
        sf: np.ndarray,
        eref=None,
        figsize=(4, 3),
        ylim=(-3, 3),
        dpi=150,
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
    ):
        """
        Plot the spectral function.

        :param np.ndarray engs: The energies of the spectral functions.
        :param sf: An array of the spectral function.
        :param float eref: Reference energy to be used - this energy will be set as zero.
        :param figsize: Size of the figure.
        :param ylim: Plotting limit for the y-axis, with respect to `eref`.
        :param dpi: DPI of the generated graph.
        :param vscale: A scaling factor for the colour map.
        :param contour_plot): Whether to use contour plot instead of normal meshed color map.
        :param alphathe color map.
        :param savethe file where the generated figure is saved.
        :param ax: An existing plotting axis to be be used.
        :param vminfor the color map.
        :param vmaxfor the color map.
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

        if ax is None:
            if nspin == 1:
                fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
                axes = [ax]
            else:
                fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
        else:
            if not isinstance(ax, list):
                axes = [ax]
            else:
                axes = ax
            fig = axes[0].figure

        # Shift the kdist so the pcolormesh draw the pixel centred on the original point
        X, Y = np.meshgrid(kdist, engs - eref)

        # Calculate the min and max values within the field of view, scaled by the factor
        mask = (engs < (ylim[1] + eref)) & (engs > (ylim[0] + eref))
        vmin = sf[:, mask, :].min() if vmin is None else vmin

        if vmax is None:
            vmax = sf[:, mask, :].max()
            vmax = (vmax - vmin) * vscale + vmin

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

        fig.tight_layout(pad=0.2)
        if save:
            fig.savefig(save, dpi=dpi)
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
        intensity: float = 1.0,
        save: bool = False,
        ax: Union[None, plt.Axes] = None,
        show: bool = False,
        title: Union[None, str] = None,
        vmin: Union[None, float] = None,
        vmax: Union[None, float] = None,
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
        :param ax: An existing plotting axis to be be used.
        :param vmin: Lower bound for constructing the alpha channel.
        :param vmax: Upper bound for constructing the alpha channel.
        :param intensity: Scaling factor for the alpha channel.
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
            if not isinstance(ax, list):
                axes = [ax]
            else:
                axes = ax
            fig = axes[0].figure

        mask = (engs < (ylim[1] + eref)) & (engs > (ylim[0] + eref))
        vmin = sf[:, mask, :, 3].min() if vmin is None else vmin

        if vmax is None:
            vmax = sf[:, mask, :, 3].max()

        # Clip the alpha range
        alpha = sf[:, :, :, 3]
        alpha = (alpha - vmin) / (vmax - vmin) * intensity
        alpha = np.clip(alpha, 0, 1)
        sf[:, :, :, 3] = alpha

        for ispin, ax_ in zip(range(nspin), axes):
            # plt.imshow's pixel coordinates are not at the center of the pixel
            # Hence, the `extent`` need to update so teh center of the pixel is aligned with the coordinates.
            extent = np.array([0, sf.shape[2], max(engs) - eref, min(engs) - eref])
            ebin = (max(engs) - min(engs)) / sf.shape[1]
            extent[:2] -= 0.5
            extent[2:] -= ebin * 0.5
            ax_.imshow(sf[ispin], extent=extent, aspect='auto', origin='upper')
            ax_.set_ylim(ylim)
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
        labels = self.unfold.kpoint_labels
        kdist = self.unfold.get_kpoint_distances()

        # Explicit label indices
        tick_locs = []
        tick_labels = []
        for index, label in labels:
            if x_is_kidx:
                xloc = index
            else:
                xloc = kdist[index]
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
        xwidth = 0.2
        if eref is None:
            eref = self.unfold.calculated_quantities['vbm']

        fig, axes = plt.subplots(1, len(all_k))
        if effective_mass_data is None:
            effective_mass_data = eff.get_effective_masses()

        # Plot the spectral function
        for (ik, ax) in zip(all_k, axes):
            self.plot_spectral_function(engs, sf, ax=ax, eref=eref, **kwargs)
            xk = kdist[ik]
            xlim = (xk - xwidth / 2, xk + xwidth / 2)
            ax.set_xlim(xlim)
            ax.set_title(f'Kpoint: {ik}')

        elec = effective_mass_data['electrons']
        # Plot the detected effective mass fitting data on top
        for entry in elec:
            ik = entry['kpoint_index']
            iax = all_k.index(ik)
            x = entry['raw_data']['kpoint_distances']
            y = entry['raw_data']['effective_energies']
            axes[iax].plot(x, np.asarray(y) - eref, '-o', color='C1')

        hole = effective_mass_data['holes']
        for entry in hole:
            ik = entry['kpoint_index']
            iax = all_k.index(ik)
            x = entry['raw_data']['kpoint_distances']
            y = entry['raw_data']['effective_energies']
            axes[iax].plot(x, np.asarray(y) - eref, '-o', color='C2')
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
        :param ax: Existing plotting axes to be be used (list if having two spin channels).
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
            if not isinstance(ax, list):
                axes = [ax]
            else:
                axes = ax
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
            procar: Union[str, list],
            eref=None,
            gamma=False,
            npoints=2000,
            sigma=0.2,
            ncl=False,
            symm_average=True,
            figsize=(4, 3),
            ylim=(-3, 3),
            dpi=150,
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
            atoms_idx=None,
            orbitals=None,
            intensity=1.0,
            use_subplot=False,
            colors=('r', 'g', 'b', 'purple'),
            colorspace='lab',
    ):
        """
        Plot projected spectral function onto multiple subplots or a single plot with color mapping.

        This simply computes the spectral function at each orbital/atoms sites and plot them onto
        multiple subplots. The columns are for each orbital and the rows are for each spin channel.

        :param procar: Name of names of the `PROCAR` files.

        :returns: Generated plot.
        """

        unfoldset = self.unfold
        unfoldset.load_procar(procar)
        nspins = unfoldset.calculated_quantities['spectral_weights_per_set'][0].shape[0]

        atoms_idx_subplots = atoms_idx.split('|')
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
            # Setup the atoms_idx and orbitals
            this_idx, this_orbitals = process_projection_options(this_idx, this_orbitals)
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
            vmax = (vmax - vmin) * vscale + vmin
            vmaxs.append(vmax)

        # Workout the vmax and vmin
        vmax = max(vmaxs)
        vmin = 0.

        if use_subplot:
            fig, axs = plt.subplots(nspins, len(atoms_idx_subplots), sharex=True, sharey=True, squeeze=False, figsize=(3.0 * nsub, 4.0))
            # Plot the spectral function with constant colour scales
            for spectral_function, i in zip(all_sf, range(len(atoms_idx_subplots))):
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
                )
        else:
            # Make a combined plot
            sf_size = all_sf[0].shape
            stacked_sf = np.stack(all_sf, axis=-1).reshape(np.prod(sf_size), len(all_sf))

            # Construct the color basis
            colors = colors[:len(all_sf)]
            # Compute spectral weight data with RGB reshape it back into the shape (nengs, nk, 3)
            sf_rgb = interpolate_colors(colors, stacked_sf, colorspace, normalize=True).reshape(sf_size + (3,))
            sf_sum = np.sum(all_sf, axis=0)[:, :, :, None]
            sf_rgba = np.concatenate([sf_rgb, sf_sum], axis=-1)

            if ax is None:
                fig, axs = plt.subplots(1, nspins, figsize=figsize, squeeze=True)

            self._plot_spectral_function_rgba(
                eng,
                sf_rgba,
                eref=eref,
                save=save,
                title=title,
                intensity=intensity,
                ax=axs,
                figsize=figsize,
                show=show,
                dpi=dpi,
                ylim=ylim,
            )

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

        :return: A figure with plotted data.
        """
        data = efm.get_effective_masses(npoints=npoints)[carrier][idx]

        x = data['raw_data']['kpoint_distances']
        y = data['raw_data']['effective_energies']
        me = data['effective_mass']
        x1, y1 = fitted_band(x, me)
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.plot(x, y, 'x-', label='Energy ')
        ax.plot(x1, y1 + y[0], label='fitted')
        ax.legend()
        if save:
            fig.savefig(save, dpi=dpi)
        return fig


def interpolate_colors(colors: Sequence, weights: list, colorspace='lab', normalize=True):
    """
    Interpolate colors at a number of points within a colorspace.

    :param colors: A list of colors specified in any way supported by matplotlib.
    :param weights: A list of weights with the shape (n, N). Where the N values of
      the last axis give the amount of N colors supplied in `colors`.
    :param colorspace: The colorspace in which to perform the interpolation. The
            allowed values are rgb, hsv, lab, luvlc, lablch, and xyz.

    :returns: A list of colors, specified in the rgb format as a (n, 3) array.
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
    colors_srgb = [sRGBColor(*to_rgb(c)) for c in colors]

    colors_basis = [np.array(convert_color(srgb, colorspace, target_illuminant='d50').get_value_tuple()) for srgb in colors_srgb]

    # ensure weights is a numpy array
    weights = np.asarray(weights)

    # Normalise the weights if needed
    if normalize:
        weights = weights / np.linalg.norm(weights, axis=1)[:, None]

    # perform the interpolation in the colorspace basis
    interpolated_colors = colors_basis[0] * weights[:, 0][:, None]
    for i in range(1, len(colors_basis)):
        interpolated_colors += colors_basis[i] * weights[:, i][:, None]

    # convert the interpolated colors back to RGB
    rgb_colors = [convert_color(colorspace(*c), sRGBColor).get_value_tuple() for c in interpolated_colors]
    rgb_colors = np.stack(rgb_colors, axis=0)

    # ensure all rgb values are less than 1 (sometimes issues in interpolation gives
    np.clip(rgb_colors, 0, 1, rgb_colors)
    return rgb_colors
