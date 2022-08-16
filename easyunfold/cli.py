"""
Commandline interface
"""
from pathlib import Path
from monty.serialization import loadfn
import numpy as np
import click
from ase.io import read

from easyunfold.unfold import create_white_colormap_from_existing, parse_atoms_idx

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
    click.option('--npoints', type=int, default=2000, help='Number of bins for the energy.')(func)
    click.option('--sigma', type=float, default=0.02, help='Smearing width for the energy in eV.')(func)
    click.option('--eref', type=float, help='Reference energy in eV.')(func)
    click.option('--emin', type=float, help='Minimum energy in eV relative to the reference.')(func)
    click.option('--emax', type=float, help='Maximum energy in eV relative to the reference.')(func)
    click.option('--vscale', type=float, help='A scaling factor for the colour mapping.', default=1.0)(func)
    click.option('--out-file', default='unfold.png', help='Name of the output file.')(func)
    click.option('--cmap', default='PuRd', help='Name of the colour map to use.')(func)
    click.option('--show', is_flag=True, default=False, help='Show the plot interactively.')(func)
    click.option('--no-symm-average', is_flag=True, default=False, help='Do not include symmetry related kpoints for averaging.')(func)
    click.option('--procar', multiple=True, help='PROCAR files used for atomic weighting')(func)
    click.option('--atoms-idx', help='Indices of the atoms to be used for weighting (1-indexed).')(func)
    click.option('--orbitals', help='Orbitals of to be used for weighting.')(func)
    click.option('--title', help='Title to be used')(func)
    return func


@unfold.command('effective-mass')
@click.pass_context
@click.option('--intensity-tol', type=float, default=0.1)
@click.option('--spin', type=int, default=0)
@click.option('--npoints', type=int, default=3)
@click.option('--edge-detect-tol', type=float, default=0.1)
@click.option('--nocc', type=int)
def unfold_effective_mass(ctx, intensity_tol, spin, npoints, edge_detect_tol, nocc):
    """
    Compute and print effective masses by tracing the unfolded weights.

    Note that this functionality only works for simple unfolded band structures,
    and it is likely to fail for complex cases.
    """
    from easyunfold.effective_mass import EffectiveMass
    from easyunfold.unfold import UnfoldKSet
    unfoldset: UnfoldKSet = ctx.obj['obj']
    efm = EffectiveMass(unfoldset, intensity_tol=intensity_tol, edge_detect_tol=edge_detect_tol)
    click.echo(f'Band extrema data: {efm.get_band_extrema()}')
    if nocc:
        efm.set_nocc(nocc)
    output = efm.get_effective_masses(ispin=spin, npoints=npoints)

    ## Print data
    def print_data(entries, tag='me'):
        """Print the effective mass data"""
        for entry in entries:
            me = entry['effective_mass']
            kf = entry['kpoint_from']
            lf = entry['kpoint_label_from']
            kt = entry['kpoint_to']
            lt = entry['kpoint_label_to']
            click.echo(f'   {tag}: {me:.3f} {kf} ({lf}) -> {kt} ({lt})')

    click.echo('Electron effective masses:')
    print_data(output['electrons'], 'm_e')
    print('')
    click.echo('Hole effective masses:')
    print_data(output['holes'], 'm_h')

    click.echo('Unfolded band structure can be ambiguous, please cross-check with the spectral function plot.')


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
@click.option('--atoms-idx', help='Indices of the atoms to be used for weighting', required=True)
@click.option('--procar', multiple=True, help='PROCAR files used for atomic weighting', required=True)
@click.option('--combined/--no-combined', is_flag=True, default=False, help='Plot all projections in a combined graph.')
@click.option('--cmap', default='PuRd', help='Name of the colour map(s) to use. Passing a list separated by "|" for the combined plot.')
def unfold_plot_projections(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale, procar,
                            atoms_idx, orbitals, title, combined):
    """
    Plot with subplots with multiple atomic projections
    """
    from easyunfold.unfold import UnfoldKSet, EBS_cmaps
    import matplotlib.pyplot as plt

    unfoldset: UnfoldKSet = ctx.obj['obj']
    click.echo(f'Loading projections from: {procar}')
    unfoldset.load_procar(procar)

    atoms_idx_subplots = atoms_idx.split('|')
    if orbitals is not None:
        orbitals_subsplots = atoms_idx.split('|')

        # Special case: if only one set is passed, apply it to all atomic specifications
        if len(orbitals_subsplots) == 1:
            orbitals_subsplots = orbitals_subsplots * len(atoms_idx_subplots)

        if len(orbitals_subsplots) != len(atoms_idx_subplots):
            click.echo('Please pass same number of orbital specifications as the atomic indices')
            raise click.Abort()
    # If not set, use all for all subsets
    else:
        orbitals_subsplots = ['all'] * len(atoms_idx_subplots)

    # Load the data

    nsub = len(atoms_idx_subplots)

    all_sf = []
    vmaxs = []

    if eref is None:
        eref = unfoldset.calculated_quantities.get('vbm', 0.0)
    click.echo(f'Using a reference energy of {eref:.3f} eV')

    # Collect spectral functions and scale
    for this_idx, this_orbitals in zip(atoms_idx_subplots, orbitals_subsplots):
        # Setup the atoms_idx and orbitals
        this_idx, this_orbitals = process_projection_options(this_idx, this_orbitals)
        eng, spectral_function = unfoldset.get_spectral_function(gamma=gamma,
                                                                 npoints=npoints,
                                                                 sigma=sigma,
                                                                 ncl=ncl,
                                                                 atoms_idx=this_idx,
                                                                 orbitals=this_orbitals,
                                                                 symm_average=not no_symm_average)
        all_sf.append(spectral_function)

        if emin is None:
            emin = eng.min() - eref
        if emax is None:
            emax = eng.max() - eref

        # Clip the effective range
        mask = (eng < (emax + eref)) & (eng > (emin + eref))
        vmin = spectral_function[:, mask, :].min()
        vmax = spectral_function[:, mask, :].max()
        vmax = (vmax - vmin) * vscale + vmin
        vmaxs.append(vmax)

    # Workout the vmax and vmin
    vmax = max(vmaxs)
    vmin = 0.

    if not combined:
        fig, axs = plt.subplots(1, len(atoms_idx_subplots), sharex=True, sharey=True, squeeze=False, figsize=(3.0 * nsub, 4.0))
        # Plot the spectral function with constant colour scales
        for spectral_function, ax in zip(all_sf, axs[0]):
            _ = EBS_cmaps(
                unfoldset.kpts_pc,
                unfoldset.pc_latt,
                eng,
                spectral_function,
                eref=eref,
                save=out_file,
                show=False,
                explicit_labels=unfoldset.kpoint_labels,
                ylim=(emin, emax),
                vscale=vscale,
                vmax=vmax,
                vmin=vmin,
                cmap=cmap,
                title=title,
                ax=ax,
            )
    else:
        fig, ax = plt.subplots(1, 1, figsize=(3.0, 4.0))
        sf_sum = sum(all_sf)
        if '|' in cmap:
            cmaps = cmap.split('|')
        else:
            cmaps = DEFAULT_CMAPS

        for spectral_function, cmap_ in zip(all_sf, cmaps):
            cmap_ = create_white_colormap_from_existing(cmap_)
            alpha_mask = spectral_function / sf_sum * 1.5
            alpha_mask = alpha_mask.clip(0, 1.0)
            _ = EBS_cmaps(unfoldset.kpts_pc,
                          unfoldset.pc_latt,
                          eng,
                          spectral_function,
                          eref=eref,
                          save=out_file,
                          show=False,
                          explicit_labels=unfoldset.kpoint_labels,
                          ylim=(emin, emax),
                          vscale=vscale,
                          vmax=vmax,
                          vmin=vmin,
                          cmap=cmap_,
                          title=title,
                          ax=ax,
                          alpha=alpha_mask)

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
    from easyunfold.unfold import UnfoldKSet, EBS_cmaps

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

    _ = EBS_cmaps(
        unfoldset.kpts_pc,
        unfoldset.pc_latt,
        eng,
        spectral_function,
        eref=eref,
        save=out_file,
        show=False,
        explicit_labels=unfoldset.kpoint_labels,
        ylim=(emin, emax),
        vscale=vscale,
        cmap=cmap,
        title=title,
        ax=ax,
    )
    if out_file:
        click.echo(f'Unfolded band structure saved to {out_file}')

    if show:
        import matplotlib.pyplot as plt
        plt.show()


def print_symmetry_data(kset):
    """Print the symmetry information"""

    # Print space group information
    sc_spg = kset.metadata['symmetry_dataset_sc']
    click.echo('Primitive cell information:')
    click.echo(' ' * 8 + f'Space group number: {sc_spg["number"]}')
    click.echo(' ' * 8 + f'Internation symbol: {sc_spg["international"]}')
    click.echo(' ' * 8 + f'Point group: {sc_spg["pointgroup"]}')

    pc_spg = kset.metadata['symmetry_dataset_pc']
    click.echo('\nSupercell cell information:')
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


def process_projection_options(atoms_idx, orbitals):
    """Process commandline type specifications"""
    indices = parse_atoms_idx(atoms_idx)
    if orbitals and orbitals != 'all':
        orbitals = [token.strip() for token in orbitals.split(',')]
    else:
        orbitals = 'all'
    return indices, orbitals
