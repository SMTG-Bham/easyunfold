"""
Plotting utilities
"""
import matplotlib.pyplot as plt
import numpy as np
from .unfold import UnfoldKSet, clean_latex_string
from .effective_mass import EffectiveMass

# pylint: disable=too-many-locals


class UnfoldPlotter:
    """A collection plotting tools for unfolded band structures"""

    def __init__(self, unfold: UnfoldKSet):
        """Instantiate the plotter object"""
        self.unfold = unfold

    def plot_spectral_function(
        self,
        engs,
        sf,
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
        Plot spectral function.

        Args:
            engs (numpy.ndarray): The energies of the spectral functions.
            sf (np.ndarray): An array of the spectral function.
            eref (float): Reference energy to be used - this energy will be set as zero.
            figsize: Size of the figure.
            ylim: Plotting limit for the y-axis, with respect to `eref`.
            dpi: DPI of the generated graph.
            vscale: A scaling factor for the colour map.
            contour_plot (bool): Whether to use contour plot instead of normal meshed color map.
            alpha (float): Alpha for the color map.
            save (str): Name of the file where the generated figure is saved.
            ax: An existing plotting axis to be be used.
            vmin (float): Minimum value for the color map.
            vmax (float): Maximum value for the color map.
            cmap (str): Name of the color map to be used.

        Returns:
            The figure generated containing the spectral function.
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
            fig.savefig(save, dpi=300)
        if show:
            fig.show()
        return fig

    def _add_kpoint_labels(self, ax):
        """Add labels to the kpoints"""
        # Label the kpoints
        labels = self.unfold.kpoint_labels
        kdist = self.unfold.get_kpoint_distances()

        # Explicit label indices
        tick_locs = []
        tick_labels = []
        for index, label in labels:
            xloc = kdist[index]
            ax.axvline(x=xloc, lw=0.5, color='k', ls=':', alpha=0.8)
            tick_locs.append(xloc)
            tick_labels.append(clean_latex_string(label))
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labels)

    def plot_effective_mass(self, eff: EffectiveMass, engs, sf, eref=None, **kwargs):
        """
        Plot the effective masses on top of the spectral function.

        Args:
            eff: An `EffectiveMass` object used for plotting.
            engs (numpy.ndarray): The energies of the spectral functions.
            sf (np.ndarray): An array of the spectral function.
            eref (float): Reference energy to be used - this energy will be set as zero.
            **kwargs: Other keyword arguments supplied to `plot_spectral_function`.

        Returns:
            A figure with the data used for fitting effective mass plotted on top of the spectral function.
        """

        kcbm = eff.get_band_extrema(mode='cbm')[0]
        kvbm = eff.get_band_extrema(mode='vbm')[0]
        all_k = sorted(list(set(list(kcbm) + list(kvbm))))
        kdist = self.unfold.get_kpoint_distances()
        xwidth = 0.2
        if eref is None:
            eref = self.unfold.calculated_quantities['vbm']

        fig, axes = plt.subplots(1, len(all_k))
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

        return fig

    def plot_spectral_weights(
            self,
            figsize=(4, 3),
            ylim=(-3, 3),
            dpi=150,
            factor=3.0,
            eref=None,
            color='C1',
            alpha=0.5,
            save=None,
            ax=None,
            show=False,
            title=None,
    ):
        """
        Plot the spectral weights.

        Note: The reduction of symmetry means there can be multiple supercell kpoints for each
        primitive cell kpoint. When using this scattering plot representation, the markers can
        overlap each other leading to misrepresentations of the actual effective band structure.

        However, this function is still useful when: 1. the symmetry splitting is turned off. 2.
        direct visualization of the underlying spectral weight is needed. 3. Check the correctness
        of effective mass extraction.

        Args:
            eref (float): Reference energy to be used - this energy will be set as zero.
            figsize: Size of the figure.
            ylim: Plotting limit for the y-axis, with respect to `eref`.
            dpi: DPI of the generated graph.
            alpha (float): Alpha for the markers.
            save (str): Name of the file where the generated figure is saved.
            ax: Existing plotting axes to be be used (list if having two spin channels).
            factor (float): Scaling factor for the marker size.
            color (str): Color of the markers.

        Returns:
            A Figure with the spectral weights plotted as a scatter plot.
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
            fig.savefig(save, dpi=300)
        if show:
            fig.show()
        return fig
