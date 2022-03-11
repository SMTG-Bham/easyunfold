"""
Commandline interface
"""
from pathlib import Path
from monty.serialization import loadfn
import numpy as np
import click
from ase.io import read

# pylint:disable=import-outside-toplevel


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
@click.option('--nk-per-split', help='Number of kpoints per split.', type=int)
@click.option('--scf-kpoints',
              help='File (IBZKPT) to provide SCF kpoints for self-consistent calculations.',
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
    print_symmetry_data(unfoldset)

    out_file = Path(out_file)
    if scf_kpoints is not None:
        scf_kpt, _, _, scf_weights = read_kpoints(scf_kpoints)
        unfoldset.write_sc_kpoints(f'KPOINTS_{out_file.stem}', nk_per_split=nk_per_split, scf_kpoints_and_weights=(scf_kpt, scf_weights))
    else:
        unfoldset.write_sc_kpoints(f'KPOINTS_{out_file.stem}', nk_per_split=nk_per_split)
    click.echo('Supercell kpoints written to KPOITNS_' + out_file.stem)

    # Serialize the data
    Path(out_file).write_text(unfoldset.to_json())

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

    click.echo(f'\nNo. of k points in the primitive cell         : {unfoldset.nkpts_orig}')
    click.echo(f'No. of expanded kpoints to be calculated cell   : {unfoldset.nkpts_expand}')
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
    """Perform the unfolding"""
    from easyunfold.unfold import UnfoldKSet

    unfoldset: UnfoldKSet = ctx.obj['obj']
    unfoldset.get_spectral_weights(wavecar, gamma, ncl=ncl)

    out_path = save_as if save_as else ctx.obj['fname']
    Path(out_path).write_text(unfoldset.to_json())
    click.echo('Unfolding data written to ' + out_path)


@unfold.command('plot')
@click.pass_context
@click.option('--gamma', is_flag=True, help='Is the calculation a gamma only one?')
@click.option('--ncl', is_flag=True, help='Is the calculation with non-colinear spin?')
@click.option('--npoints', type=int, default=2000, help='Number of bins for the energy.')
@click.option('--sigma', type=float, default=0.02, help='Smearing width for the energy in eV.')
@click.option('--eref', type=float, help='Reference energy in eV.')
@click.option('--emin', type=float, help='Minimum energy in eV relative to the reference.')
@click.option('--emax', type=float, help='Maximum energy in eV relative to the reference.')
@click.option('--vscale', type=float, help='A scaling factor for the colour mapping.', default=1.0)
@click.option('--out-file', default='unfold.png', help='Name of the output file.')
@click.option('--cmap', default='PuRd', help='Name of the colour map to use.')
@click.option('--show', is_flag=True, default=False, help='Show the plot interactively.')
@click.option('--no-symm-average', is_flag=True, default=False, help='Do not include symmetry related kpoints for averaging.')
def unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl, no_symm_average, vscale):
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

    eng, spectral_function = unfoldset.get_spectral_function(gamma=gamma,
                                                             npoints=npoints,
                                                             sigma=sigma,
                                                             ncl=ncl,
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
