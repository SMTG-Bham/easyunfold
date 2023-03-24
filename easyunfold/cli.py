"""
Commandline interface
"""
from pathlib import Path
from monty.serialization import loadfn
import numpy as np
import click
from ase.io import read

from easyunfold.unfold import process_projection_options

# pylint:disable=import-outside-toplevel,too-many-locals

DEFAULT_CMAPS = [
    'Purples', 'Greens', 'Oranges', 'Reds', 'Blue', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn',
    'BuGn', 'YlGn'
]


@click.group('easyunfold')
def easyunfold():
    """
    Tool for performing band unfolding
    """
    return


@easyunfold.command()
@click.option('--time-reversal/--no-time-reversal', default=True)
@click.argument('pc-file')
@click.argument('sc-file')
@click.argument('kpoints')
@click.option('--matrix', '-m', help='Transformation matrix')
@click.option('--symprec', help='Transformation matrix', type=float, default=1e-5)
@click.option('--out-file', default='easyunfold.json', help='Name of the output file')
@click.option('--no-expand', help='Do not expand the kpoints by symmetry', default=False, is_flag=True)
@click.option('--nk-per-split', help='Number of band structure kpoints per split.', type=int)
@click.option('--scf-kpoints',
              help='File (IBZKPT) to provide SCF kpoints for self-consistent calculations. Needed for hybrid functional calculations.',
              type=click.Path(exists=True, dir_okay=False))
def generate(pc_file, sc_file, matrix, kpoints, time_reversal, out_file, no_expand, symprec, nk_per_split, scf_kpoints):
    """
    Generate the kpoints for sampling the supercell
    """
    from easyunfold.unfold import UnfoldKSet, read_kpoints

    primitive = read(pc_file)
    supercell = read(sc_file)
    if matrix:
        transform_matrix = matrix_from_string(matrix)
        if not np.allclose(primitive.cell @ transform_matrix, supercell.cell):
            click.echo('Warning: the super cell and the the primitive cell are not commensure.')
            click.echo('Proceed with the assumed tranform matrix')
        click.echo(f'Transform matrix:\n{transform_matrix.tolist()}')
    else:
        tmp = supercell.cell @ np.linalg.inv(primitive.cell)
        transform_matrix = np.rint(tmp)
        if not np.allclose(tmp, transform_matrix):
            click.echo('The super cell and the the primitive cell are not commensure.')
            raise click.Abort()

        click.echo(f'(Guessed) Transform matrix:\n{transform_matrix.tolist()}')

    kpoints, _, labels, _ = read_kpoints(kpoints)
    click.echo(f'{len(kpoints)} kpoints specified along the path')

    unfoldset = UnfoldKSet.from_atoms(transform_matrix,
                                      kpoints,
                                      primitive,
                                      supercell,
                                      time_reversal=time_reversal,
                                      expand=not no_expand,
                                      symprec=symprec)
    unfoldset.kpoint_labels = labels
    try:
        print_symmetry_data(unfoldset)
    except KeyError:
        pass

    out_file = Path(out_file)
    if scf_kpoints is not None:
        scf_kpt, _, _, scf_weights = read_kpoints(scf_kpoints)
        unfoldset.write_sc_kpoints(f'KPOINTS_{out_file.stem}', nk_per_split=nk_per_split, scf_kpoints_and_weights=(scf_kpt, scf_weights))
    else:
        unfoldset.write_sc_kpoints(f'KPOINTS_{out_file.stem}', nk_per_split=nk_per_split)
    click.echo('Supercell kpoints written to KPOITNS_' + out_file.stem)

    # Serialize the data
    Path(out_file).write_text(unfoldset.to_json(), encoding='utf-8')

    click.echo('Unfolding settings written to ' + str(out_file))


@easyunfold.group('unfold')
@click.option('--data-file', default='easyunfold.json', type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.pass_context
def unfold(ctx, data_file):
    """Perform unfolding and plotting"""

    unfoldset = loadfn(data_file)
    click.echo(f'Loaded data from {data_file}')
    ctx.obj = {'obj': unfoldset, 'fname': data_file}


@unfold.command('status')
@click.pass_context
def unfold_status(ctx):
    """Print the status"""
    from easyunfold.unfold import UnfoldKSet
    unfoldset: UnfoldKSet = ctx.obj['obj']
    print_symmetry_data(unfoldset)
    nkpts_sc = len(unfoldset.expansion_results['reduced_sckpts'])
    click.echo()
    click.echo(f'No. of k points in the primitive cell           : {unfoldset.nkpts_orig}')
    click.echo(f'No. of expanded kpoints to be calculated cell   : {nkpts_sc} ({unfoldset.nkpts_expand})')
    click.echo(f'No. of rotations in the primitive cell          : {unfoldset.pc_opts.shape[0]}')
    click.echo(f'No. of rotations in the super cell              : {unfoldset.sc_opts.shape[0]}')
    click.echo()
    click.echo('Path in the primitive cell:')
    for index, label in unfoldset.kpoint_labels:
        click.echo(f'   {label:<10}: {index+1:<5}')

    if unfoldset.is_calculated:
        click.echo('Unfolding had been performed - use `unfold plot` to plot the spectral function.')
    else:
        click.echo('Please run the supercell band structure calculation and run `unfold calculate`.')


@unfold.command('calculate')
@click.pass_context
@click.argument('wavecar', type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.option('--save-as')
@click.option('--gamma', is_flag=True)
@click.option('--ncl', is_flag=True)
def unfold_calculate(ctx, wavecar, save_as, gamma, ncl):
    """
    Perform the unfolding

    Multiple WAVECAR files can be supplied for split-path calculations.
    """
    from easyunfold.unfold import UnfoldKSet

    unfoldset: UnfoldKSet = ctx.obj['obj']
    unfoldset.get_spectral_weights(wavecar, gamma, ncl=ncl)

    out_path = save_as if save_as else ctx.obj['fname']
    Path(out_path).write_text(unfoldset.to_json(), encoding='utf-8')
    click.echo('Unfolding data written to ' + out_path)


def add_plot_options(func):
    """Added common plotting options to a function"""
    click.option('--gamma', is_flag=True, help='Is the calculation a gamma only one?')(func)
    click.option('--ncl', is_flag=True, help='Is the calculation with non-colinear spin?')(func)
    click.option('--npoints', type=int, default=2000, help='Number of bins for the energy.', show_default=True)(func)
    click.option('--sigma', type=float, default=0.02, help='Smearing width for the energy in ' 'eV.', show_default=True)(func)
    click.option('--eref', type=float, help='Reference energy in eV.')(func)
    click.option('--emin', type=float, default=-5., help='Minimum energy in eV relative to the ' 'reference.', show_default=True)(func)
    click.option('--emax', type=float, default=5., help='Maximum energy in eV relative to the ' 'reference.', show_default=True)(func)
    click.option('--vscale', type=float, help='A scaling factor for the colour mapping.', default=1.0)(func)
    click.option('--out-file', default='unfold.png', help='Name of the output file.')(func)
    click.option('--cmap',
                 default='PuRd',
                 help='Name of the colour map(s) to use. Passing a list separated by "|" for the '
                 'combined plot.',
                 show_default=True)(func)
    click.option('--show', is_flag=True, default=False, help='Show the plot interactively.')(func)
    click.option('--no-symm-average',
                 is_flag=True,
                 default=False,
                 help='Do not include symmetry '
                 'related kpoints for '
                 'averaging.',
                 show_default=True)(func)
    click.option('--procar', multiple=True, help='PROCAR files used for atomic weighting')(func)
    click.option('--atoms-idx', help='Indices of the atoms to be used for weighting (1-indexed).')(func)
    click.option('--orbitals', help='Orbitals of to be used for weighting.')(func)
    click.option('--title', help='Title to be used')(func)
    return func


@unfold.command('effective-mass')
@click.pass_context
@click.option('--intensity-threshold', type=float, default=0.1, help='Intensity threshold for detecting valid bands.')
@click.option('--spin', type=int, default=0, help='Index of the spin channel.')
@click.option('--npoints', type=int, default=3, help='Number of kpoints used for fitting from the extrema.')
@click.option('--extrema-detect-tol', type=float, default=0.01, help='Tolerance for band extrema detection.')
@click.option('--degeneracy-detect-tol', type=float, default=0.01, help='Tolerance for band degeneracy detection at extrema.')
@click.option('--nocc', type=int)
@click.option('--plot', is_flag=True, default=False)
@click.option('--plot-fit', is_flag=True, default=False, help='Generate plots of the band edge and parabolic fits.')
@click.option('--fit-label', help='Which branch to use for plot fitting. e.g. electrons:0', default='electrons:0')
@click.option('--band-filter', default=None, type=int, help='Only displace information for this band.')
@click.option('--out-file', default='unfold-effective-mass.png', help='Name of the output file.')
def unfold_effective_mass(ctx, intensity_threshold, spin, band_filter, npoints, extrema_detect_tol, degeneracy_detect_tol, nocc, plot,
                          plot_fit, fit_label, out_file):
    """
    Compute and print effective masses by tracing the unfolded weights.

    Note that this functionality only works for simple unfolded band structures,
    and it is likely to fail for complex cases.

    Use the `--plot-fit` flag to generate the detected band edge verses the parabolic fit
    for visual confirmation if unsure.
    """
    from easyunfold.effective_mass import EffectiveMass
    from easyunfold.unfold import UnfoldKSet
    from tabulate import tabulate
    unfoldset: UnfoldKSet = ctx.obj['obj']
    efm = EffectiveMass(unfoldset, intensity_tol=intensity_threshold, extrema_tol=extrema_detect_tol, degeneracy_tol=degeneracy_detect_tol)

    click.echo('Band extrema data:')
    table = []
    for mode in ['cbm', 'vbm']:
        for kid, subkid, iband in zip(*efm.get_band_extrema(mode=mode)):
            band_idx = ','.join(map(str, iband))
            table.append([kid, mode, subkid, band_idx])
    click.echo(tabulate(table, headers=['Kpoint index', 'Kind', 'Sub-kpoint index', 'Band indices']))
    click.echo('')

    if nocc:
        efm.set_nocc(nocc)
    output = efm.get_effective_masses(ispin=spin, npoints=npoints)

    # Filter by band if requested
    if band_filter is not None:
        for carrier in ['electrons', 'holes']:
            output[carrier] = [entry for entry in output[carrier] if entry['band_index'] == band_filter]

    ## Print data
    def print_data(entries, tag='me'):
        """Print the effective mass data"""
        table = []
        for i, entry in enumerate(entries):
            me = entry['effective_mass']

            kf = entry['kpoint_from']
            lf = entry['kpoint_label_from']
            kfrom = f'{kf} ({lf})'

            kt = [round(x, 5) for x in entry['kpoint_to']]
            lt = entry['kpoint_label_to']
            kto = f'{kt} ({lt})'
            bidx = entry['band_index']
            table.append([i, tag, me, bidx, kfrom, kto])
        click.echo(tabulate(table, headers=['index', 'Kind', 'Effective mass', 'Band index', 'from', 'to']))

    click.echo('Electron effective masses:')
    print_data(output['electrons'], 'm_e')
    print('')
    click.echo('Hole effective masses:')
    print_data(output['holes'], 'm_h')

    click.echo('Unfolded band structure can be ambiguous, please cross-check with the spectral function plot.')

    if plot:
        from easyunfold.plotting import UnfoldPlotter
        plotter = UnfoldPlotter(unfoldset)
        click.echo('Generating spectral function plot for visualising detected branches...')
        engs, sf = unfoldset.get_spectral_function()
        plotter.plot_effective_mass(efm, engs, sf, effective_mass_data=output, save=out_file)

    elif plot_fit:
        from easyunfold.plotting import UnfoldPlotter
        plotter = UnfoldPlotter(unfoldset)
        carrier, idx = fit_label.split(':')
        fname = Path(out_file).stem
        ext = Path(out_file).suffix
        for carrier in ['electrons', 'holes']:
            for idx, _ in enumerate(output[carrier]):
                plotter.plot_effective_mass_fit(efm=efm,
                                                npoints=npoints,
                                                carrier=carrier,
                                                idx=int(idx),
                                                save=fname + f'_fit_{carrier}_{idx}' + ext)


@unfold.command('plot')
@click.pass_context
@add_plot_options
def unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar, atoms_idx,
                orbitals, title):
    """
    Plot the spectral function

    This command uses the stored unfolding data to plot the effective bands structure (EBS).
    """
    _unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar, atoms_idx,
                 orbitals, title)


@unfold.command('plot-projections')
@click.pass_context
@add_plot_options
@click.option('--combined/--no-combined', is_flag=True, default=False, help='Plot all projections in a combined graph.')
@click.option('--intensity', default=1.0, help='Color intensity for combined plot', type=float, show_default=True)
@click.option('--colors', help='Colors to be used for combined plot. Comma separated. Default = "r,g,b".')
def unfold_plot_projections(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar,
                            atoms_idx, orbitals, title, combined, intensity, colors):
    """
    Plot with subplots with multiple atomic projections
    """
    from easyunfold.unfold import UnfoldKSet
    from easyunfold.plotting import UnfoldPlotter
    import matplotlib.pyplot as plt

    unfoldset: UnfoldKSet = ctx.obj['obj']
    click.echo(f'Loading projections from: {procar}')

    plotter = UnfoldPlotter(unfoldset)
    fig = plotter.plot_projected(procar,
                                 gamma=gamma,
                                 npoints=npoints,
                                 sigma=sigma,
                                 eref=eref,
                                 ylim=(emin, emax),
                                 cmap=cmap,
                                 ncl=ncl,
                                 symm_average=not no_symm_average,
                                 atoms_idx=atoms_idx,
                                 orbitals=orbitals,
                                 title=title,
                                 vscale=vscale,
                                 use_subplot=not combined,
                                 intensity=intensity,
                                 colors=colors.split(',') if colors is not None else None)

    if out_file:
        fig.savefig(out_file, dpi=300)
        click.echo(f'Unfolded band structure saved to {out_file}')

    if show:
        plt.show()


def _unfold_plot(ctx,
                 gamma,
                 npoints,
                 sigma,
                 eref,
                 out_file,
                 show,
                 emin,
                 emax,
                 cmap,
                 ncl,
                 no_symm_average,
                 vscale,
                 procar,
                 atoms_idx,
                 orbitals,
                 title,
                 ax=None):
    """
    Plot the spectral function

    This command uses the stored unfolding data to plot the effective bands structure (EBS).
    """
    from easyunfold.unfold import UnfoldKSet
    from easyunfold.plotting import UnfoldPlotter

    unfoldset: UnfoldKSet = ctx.obj['obj']
    if not unfoldset.is_calculated:
        click.echo('Unfolding has not been performed yet, please run `unfold calculate` command.')
        raise click.Abort()

    if eref is None:
        eref = unfoldset.calculated_quantities.get('vbm', 0.0)
    click.echo(f'Using a reference energy of {eref:.3f} eV')

    # Process the PROCAR
    if procar:
        click.echo(f'Loading projections from: {procar}')
        unfoldset.load_procar(procar)

    # Setup the atoms_idx and orbitals
    if atoms_idx:
        atoms_idx, orbitals = process_projection_options(atoms_idx, orbitals)
    else:
        atoms_idx = None
        orbitals = None

    eng, spectral_function = unfoldset.get_spectral_function(gamma=gamma,
                                                             npoints=npoints,
                                                             sigma=sigma,
                                                             ncl=ncl,
                                                             atoms_idx=atoms_idx,
                                                             orbitals=orbitals,
                                                             symm_average=not no_symm_average)
    if emin is None:
        emin = eng.min() - eref
    if emax is None:
        emax = eng.max() - eref

    plotter = UnfoldPlotter(unfoldset)
    fig = plotter.plot_spectral_function(eng,
                                         spectral_function,
                                         eref=eref,
                                         save=out_file,
                                         show=False,
                                         ylim=(emin, emax),
                                         vscale=vscale,
                                         cmap=cmap,
                                         title=title,
                                         ax=ax)

    if out_file:
        click.echo(f'Unfolded band structure saved to {out_file}')

    if show:
        fig.show()


def print_symmetry_data(kset):
    """Print the symmetry information"""

    # Print space group information
    sc_spg = kset.metadata['symmetry_dataset_sc']
    click.echo('Supercell cell information:')
    click.echo(' ' * 8 + f'Space group number: {sc_spg["number"]}')
    click.echo(' ' * 8 + f'Internation symbol: {sc_spg["international"]}')
    click.echo(' ' * 8 + f'Point group: {sc_spg["pointgroup"]}')

    pc_spg = kset.metadata['symmetry_dataset_pc']
    click.echo('\nPrimitive cell information:')
    click.echo(' ' * 8 + f'Space group number: {pc_spg["number"]}')
    click.echo(' ' * 8 + f'Internation symbol: {pc_spg["international"]}')
    click.echo(' ' * 8 + f'Point group: {pc_spg["pointgroup"]}')


def matrix_from_string(string):
    """Parse transform matrix from a string"""
    elems = [float(x) for x in string.split()]
    # Try gussing the transform matrix
    if len(elems) == 3:
        transform_matrix = np.diag(elems)
    else:
        transform_matrix = np.array(elems).reshape((3, 3))
    return transform_matrix
