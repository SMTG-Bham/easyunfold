"""
Commandline interface
"""
import warnings
from pathlib import Path
from monty.serialization import loadfn
import numpy as np
import click
from ase.io import read

from easyunfold.unfold import process_projection_options, parse_atoms

# pylint:disable=import-outside-toplevel,too-many-locals,too-many-arguments

SUPPORTED_DFT_CODES = ('vasp', 'castep')

DEFAULT_CMAPS = [
    'Purples', 'Greens', 'Oranges', 'Reds', 'Blue', 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn',
    'BuGn', 'YlGn'
]

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help']}

warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL fields")  # unnecessary pymatgen potcar warnings


@click.group('easyunfold', context_settings=CONTEXT_SETTINGS)
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
@click.option(
    '--code',
    '-c',
    help='Name of the DFT code that the kpoints should be formatted to.',
    default='vasp',
    type=click.Choice(SUPPORTED_DFT_CODES),
    show_default=True,
)
@click.option('--matrix', '-m', help='Transformation matrix')
@click.option('--symprec', help='Transformation matrix', type=float, default=1e-5, show_default=True)
@click.option('--out-file', default='easyunfold.json', help='Name of the output file')
@click.option('--no-expand', help='Do not expand the kpoints by symmetry', default=False, is_flag=True)
@click.option('--nk-per-split', help='Number of band structure kpoints per split.', type=int)
@click.option('--scf-kpoints',
              help='File (IBZKPT) to provide SCF kpoints for self-consistent calculations. Needed for hybrid functional calculations.',
              type=click.Path(exists=True, dir_okay=False))
@click.option('--yes', '-y', is_flag=True, default=False, help='Skip and confirmation.')
def generate(pc_file, code, sc_file, matrix, kpoints, time_reversal, out_file, no_expand, symprec, nk_per_split, scf_kpoints, yes):
    """
    Generate the kpoints for performing supercell calculations.

    There are two modes of running supercell calculations:

    1. Use the generated kpoints for unfolding for non-SCF calculations, e.g. with a fixed charged density from the SCF calculation.
    2. Include the generated kpoints in SCF calculation but set their weights to zeros.

    In both cases, the kpoints can be split into multiple calculations.

    This command defaults to go down the first route where the kpoints generated are only
    those needed for the unfolding, all kpoints and written into the same file, but the user
    is free to manually split them into multiple calculations.
    The second approach can be activated with the ``--nk-per-split`` option and in the same
    time provide the ``--scf-kpoints``. This will generate a serial of kpoints files containing
    the SCF kpoints followed by the actual kpoints needed for unfolding with zero weights.
    """
    from easyunfold.unfold import UnfoldKSet
    from easyunfold.utils import read_kpoints

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

    kpoints, _, labels, _ = read_kpoints(kpoints, code=code)
    click.echo(f'{len(kpoints)} kpoints specified along the path')

    unfoldset = UnfoldKSet.from_atoms(transform_matrix,
                                      kpoints,
                                      primitive,
                                      supercell,
                                      time_reversal=time_reversal,
                                      expand=not no_expand,
                                      dft_code=code,
                                      symprec=symprec)
    unfoldset.kpoint_labels = labels
    try:
        print_symmetry_data(unfoldset)
    except KeyError:
        pass

    out_file = Path(out_file)
    if code == 'vasp':
        out_kpt_name = f'KPOINTS_{out_file.stem}'
    elif code == 'castep':
        out_kpt_name = f'{out_file.stem}_sc_kpoints.cell'
    else:
        raise RuntimeError(f'Unknown code: {code}')

    if scf_kpoints is not None:
        scf_kpt, _, _, scf_weights = read_kpoints(scf_kpoints, code=code)
        scf_kpoints_and_weights = (scf_kpt, scf_weights)
    else:
        scf_kpoints_and_weights = None

    if Path(out_kpt_name).is_file() and not yes:
        click.confirm(f'Output file {out_kpt_name} already exists, continue?', abort=True)

    unfoldset.write_sc_kpoints(
        out_kpt_name,
        nk_per_split=nk_per_split,
        scf_kpoints_and_weights=scf_kpoints_and_weights,
        source=sc_file,
    )

    click.echo(f'Supercell kpoints written to {out_kpt_name}')

    # Serialize the data
    if Path(out_file).is_file() and not yes:
        click.confirm(f'Output file {out_file} already exists, continue?', abort=True)

    Path(out_file).write_text(unfoldset.to_json(), encoding='utf-8')

    click.echo('Unfolding settings written to ' + str(out_file))


@easyunfold.group('unfold')
@click.option('--data-file', default='easyunfold.json', type=click.Path(exists=True, file_okay=True, dir_okay=False), show_default=True)
@click.pass_context
def unfold(ctx, data_file):
    """
    Perform unfolding and plotting

    There are multiple sub-command available under this command group.
    """

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
    click.echo(f'No. of supercell k points                       : {nkpts_sc}')
    click.echo(f'No. of primitive cell symmetry operations       : {unfoldset.pc_opts.shape[0]}')
    click.echo(f'No. of supercell symmetry operations            : {unfoldset.sc_opts.shape[0]}')
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
@click.argument('wavefunc', type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.option('--save-as')
@click.option('--gamma', is_flag=True, help='If the calculation is $\\Gamma$-only')
@click.option('--ncl', is_flag=True, help='If the calculation has non-collinear spins.')
def unfold_calculate(ctx, wavefunc, save_as, gamma, ncl):
    """
    Perform the unfolding

    Multiple wave function (e.g. the WAVECAR file if using VASP) files can be supplied for
    split-path calculations.
    Once the calculation is done, the wave function data are not longer needed as the
    spectral weights will be stored in the outputs ``json`` file.
    """
    from easyunfold.unfold import UnfoldKSet

    unfoldset: UnfoldKSet = ctx.obj['obj']
    unfoldset.get_spectral_weights(wavefunc, gamma, ncl=ncl)

    out_path = save_as if save_as else ctx.obj['fname']
    Path(out_path).write_text(unfoldset.to_json(), encoding='utf-8')
    click.echo('Unfolding data written to ' + out_path)


def add_plot_options(func):
    """Added common plotting options to a function"""
    click.option('--gamma', is_flag=True, help='Is the calculation a gamma only one?', show_default=True)(func)
    click.option('--ncl', is_flag=True, help='Is the calculation with non-colinear spin?', show_default=True)(func)
    click.option('--npoints', type=int, default=2000, help='Number of bins for the energy.', show_default=True)(func)
    click.option('--sigma', type=float, default=0.02, help='Smearing width for the energy in eV.', show_default=True)(func)
    click.option('--eref', type=float, help='Reference energy in eV.')(func)
    click.option('--emin', type=float, default=-5., help='Minimum energy in eV relative to the reference.', show_default=True)(func)
    click.option('--emax', type=float, default=5., help='Maximum energy in eV relative to the reference.', show_default=True)(func)
    click.option('--colour-norm', type=float, help='A normalisation/scaling factor for the colour mapping. Smaller values will increase colour intensity.',
                 default=1.0,
                 show_default=True)(func)
    click.option('--out-file', default='unfold.png', help='Name of the output file.', show_default=True)(func)
    click.option('--cmap', default='PuRd', help='Name of the colour map to use.', show_default=True)(func)
    click.option('--show', is_flag=True, default=False, help='Show the plot interactively.')(func)
    click.option('--no-symm-average',
                 is_flag=True,
                 default=False,
                 help='Do not include symmetry related kpoints for averaging.',
                 show_default=True)(func)
    click.option('--dos', help='Path to vasprun.xml(.gz) file from which to read the density of states (DOS) information. '
                               'If set, the density of states will be plotted alongside the unfolded bandstructure.')(func)
    click.option('--dos-label', type=str, help='Axis label for DOS if included.', show_default=True)(func)
    click.option('--zero-line', is_flag=True, default=False, help='Plot horizontal line at zero energy.', show_default=True)(func)
    click.option('--dos-elements', help='Elemental orbitals to plot (e.g. "C.s.p,O") for DOS if '
                                        'included. If not specified, will be set to the orbitals of `--atoms` (if '
                                        'specified)')(func)
    click.option('--dos-orbitals', help='Orbitals to split into lm-decomposed (e.g. p -> px, py, '
                                    'pz) contributions (e.g. "S.p") for DOS if included.')(func)
    click.option('--dos-atoms', help='Atoms to include (e.g. "O.1.2.3,Ru.1.2.3") for DOS if included.')(func)
    click.option('--legend-cutoff', type=float, default=3, help='Cut-off in % of total DOS that determines if a line is given a label.', show_default=True)(func)
    click.option('--gaussian', type=float, help='Standard deviation of DOS gaussian broadening in eV.',
                 default=0.05, show_default=True)(func)
    click.option('--scale', type=float, help='Scaling factor for the DOS plot.', default=1.0)(func)
    click.option('--no-total', is_flag=True, default=False, help="Don't plot the total density of states.")(func)
    click.option('--total-only', is_flag=True, default=False, help="Only plot the total density of states.")(func)
    click.option('--procar',
                 multiple=True,
                 default=["PROCAR"],
                 help=('PROCAR file(s) for atomic weighting, can be passed multiple times if more '
                       'than one PROCAR should be used. Default is to read PROCAR in current directory'))(func)
    click.option('--atoms', help='Atoms to be used for weighting, as a comma-separated list (e.g. "Na,'
                                 'Bi,S"). The POSCAR or CONTCAR file must be present in the current '
                                 'directory for this, otherwise use `--atoms-idx`.')(func)
    click.option('--atoms-idx', help='Indices of the atoms to be used for weighting (1-indexed), '
                                     'comma-separated and "-" can be used to define ranges, in groups '
                                     'separated by "|". E.g. "1-20|21,22,23"')(func)
    click.option('--orbitals', help='Orbitals to be used for weighting, comma-separated, in groups '
                                    'separated by "|". E.g. "px,py|pz"')(func)
    click.option('--title', help='Title to be used')(func)
    click.option('--width', help='Width of the figure', type=float, default=4., show_default=True)(func)
    click.option('--height', help='Height of the figure', type=float, default=3., show_default=True)(func)
    click.option('--dpi', help='DPI for the figure when saved as raster image.', type=int, default=300, show_default=True)(func)
    click.option('--mpl-style-file',
                  type=click.Path(exists=True, file_okay=True, dir_okay=False),
                  help='Use this file to customise the matplotlib style sheet')(func)
    return func


@unfold.command('effective-mass')
@click.pass_context
@click.option('--intensity-threshold', type=float, default=0.1, help='Intensity threshold for detecting valid bands.', show_default=True)
@click.option('--spin', type=int, default=0, help='Index of the spin channel.', show_default=True)
@click.option('--npoints', type=int, default=3, help='Number of kpoints used for fitting from the extrema.', show_default=True)
@click.option('--extrema-detect-tol', type=float, default=0.01, help='Tolerance for band extrema detection.', show_default=True)
@click.option('--degeneracy-detect-tol',
              type=float,
              default=0.01,
              help='Tolerance for band degeneracy detection at extrema.',
              show_default=True)
@click.option('--nocc', type=int, help='DEV: Use this band as the extrema at all kpoints.')
@click.option('--plot', is_flag=True, default=False)
@click.option('--plot-fit', is_flag=True, default=False, help='Generate plots of the band edge and parabolic fits.')
@click.option('--fit-label', help='Which branch to use for plot fitting. e.g. electrons:0', default='electrons:0', show_default=True)
@click.option('--band-filter', default=None, type=int, help='Only displace information for this band.')
@click.option('--out-file', default='unfold-effective-mass.png', help='Name of the output file.', show_default=True)
def unfold_effective_mass(ctx, intensity_threshold, spin, band_filter, npoints, extrema_detect_tol, degeneracy_detect_tol, nocc, plot,
                          plot_fit, fit_label, out_file):
    """
    Compute and print effective masses by tracing the unfolded weights.

    Note that this functionality only works for simple unfolded band structures,
    and it is likely to fail for complex cases.

    Use the ``--plot-fit`` flag to generate the detected band edge verses the parabolic fit
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
def unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl,
                no_symm_average, colour_norm, dos, dos_label, zero_line, dos_elements, dos_orbitals,
                dos_atoms, legend_cutoff, gaussian, no_total, total_only, scale,
                procar, atoms, atoms_idx, orbitals, title, width, height, dpi, mpl_style_file):
    """
    Plot the spectral function

    This command uses the stored unfolding data to plot the effective bands structure (EBS) using the spectral function.
    """
    if mpl_style_file:
        click.echo(f'Using custom plotting style from {mpl_style_file}')
        import matplotlib.style
        matplotlib.style.use(mpl_style_file)

    _unfold_plot(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl,
                 no_symm_average, colour_norm, dos, dos_label, zero_line, dos_elements, dos_orbitals,
                 dos_atoms, legend_cutoff, gaussian, no_total, total_only, scale,
                 procar, atoms, atoms_idx, orbitals, title, width, height, dpi)


@unfold.command('plot-projections')
@click.pass_context
@add_plot_options
@click.option('--combined/--no-combined', is_flag=True, default=False, help='Plot all projections in a combined graph.')
@click.option('--colors', help='Colors to be used for combined plot, comma separated.', default='r,g,b,purple', show_default=True)
def unfold_plot_projections(ctx, gamma, npoints, sigma, eref, out_file, show, emin, emax, cmap, ncl,
                            no_symm_average, colour_norm, dos, dos_label, zero_line, dos_elements, dos_orbitals,
                            dos_atoms, legend_cutoff, gaussian, no_total, total_only, scale,
                            procar, atoms, atoms_idx, orbitals, title, combined, colors, width, height, dpi, mpl_style_file):
    """
    Plot the effective band structure with atomic projections.
    """
    from easyunfold.unfold import UnfoldKSet
    from easyunfold.plotting import UnfoldPlotter
    import matplotlib.pyplot as plt

    if mpl_style_file:
        click.echo(f'Using custom plotting style from {mpl_style_file}')
        import matplotlib.style
        matplotlib.style.use(mpl_style_file)

    unfoldset: UnfoldKSet = ctx.obj['obj']
    click.echo(f'Loading projections from: {procar}')

    plotter = UnfoldPlotter(unfoldset)
    if dos:
        from sumo.plotting.dos_plotter import SDOSPlotter
        from sumo.electronic_structure.dos import load_dos
        from sumo.cli.dosplot import _el_orb, _atoms

        dos_elements = _el_orb(dos_elements) if dos_elements is not None else None
        dos_orbitals = _el_orb(dos_orbitals) if dos_orbitals is not None else None
        dos_atoms = _atoms(dos_atoms) if dos_atoms is not None else None

        # Set dos_elements to match atoms if set and dos_elements not specified
        if atoms and not dos_elements:
            parsed_atoms, _idx, _orbitals = parse_atoms(atoms, orbitals)
            dos_elements = ",".join([f"{atom}.s.p.d.f" for atom in parsed_atoms])
            dos_elements = _el_orb(dos_elements)

        warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL fields")  # unnecessary pymatgen potcar warnings
        dos, pdos = load_dos(dos,
                             dos_elements,
                             dos_orbitals,
                             dos_atoms,
                             gaussian,
                             total_only,
                             )
        dos_plotter = SDOSPlotter(dos, pdos)
        dos_options = {
            "plot_total": not no_total,
            "legend_cutoff": legend_cutoff,
            "yscale": scale,
        }
    else:
        dos_plotter = None
        dos_options = None

    fig = plotter.plot_projected(procar,
                                 dos_plotter=dos_plotter,
                                 dos_label=dos_label,
                                 dos_options=dos_options,
                                 zero_line=zero_line,
                                 gamma=gamma,
                                 npoints=npoints,
                                 sigma=sigma,
                                 eref=eref,
                                 ylim=(emin, emax),
                                 cmap=cmap,
                                 ncl=ncl,
                                 symm_average=not no_symm_average,
                                 atoms=atoms,
                                 atoms_idx=atoms_idx,
                                 orbitals=orbitals,
                                 title=title,
                                 colour_norm=colour_norm,
                                 use_subplot=not combined,
                                 figsize=(width, height),
                                 colors=colors.split(',') if colors is not None else None)

    if out_file:
        fig.savefig(out_file, dpi=dpi)
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
                 colour_norm,
                 dos,
                 dos_label,
                 zero_line,
                 dos_elements,
                 dos_orbitals,
                 dos_atoms,
                 legend_cutoff,
                 gaussian,
                 no_total,
                 total_only,
                 scale,
                 procar,
                 atoms,
                 atoms_idx,
                 orbitals,
                 title,
                 width,
                 height,
                 dpi,
                 ax=None):
    """
    Routine for plotting the spectral function.
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

    # Setup the atoms_idx and orbitals
    if atoms or atoms_idx:
        # Process the PROCAR
        click.echo(f'Loading projections from: {procar}')
        try:
            unfoldset.load_procar(procar)
        except FileNotFoundError:
            click.echo(f'Could not find and parse the --procar file: {procar} – needed for atomic projections!')
            raise click.Abort()
        if atoms_idx:
            atoms_idx, orbitals = process_projection_options(atoms_idx, orbitals)
        else:  # parse atoms
            atoms, atoms_idx, orbitals = parse_atoms(atoms, orbitals)
            atoms_idx = [idx for sublist in atoms_idx for idx in sublist]  # flatten atoms_idx for non-projected plot

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
    if dos:
        from sumo.plotting.dos_plotter import SDOSPlotter
        from sumo.electronic_structure.dos import load_dos

        warnings.filterwarnings("ignore",
                                message="No POTCAR file with matching TITEL fields")  # unnecessary pymatgen potcar warnings
        dos, pdos = load_dos(dos,
                             dos_elements,
                             dos_orbitals,
                             dos_atoms,
                             gaussian,
                             total_only,
                             )
        dos_plotter = SDOSPlotter(dos, pdos)
        dos_options = {
            "plot_total": not no_total,
            "legend_cutoff": legend_cutoff,
            "yscale": scale,
        }
    else:
        dos_plotter = None
        dos_options = None

    fig = plotter.plot_spectral_function(eng,
                                         spectral_function,
                                         dos_plotter=dos_plotter,
                                         dos_label=dos_label,
                                         dos_options=dos_options,
                                         zero_line=zero_line,
                                         eref=eref,
                                         figsize=(width, height),
                                         save=out_file,
                                         show=False,
                                         ylim=(emin, emax),
                                         colour_norm=colour_norm,
                                         cmap=cmap,
                                         title=title,
                                         dpi=dpi,
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
    click.echo(' ' * 8 + f'International symbol: {sc_spg["international"]}')
    click.echo(' ' * 8 + f'Point group: {sc_spg["pointgroup"]}')

    pc_spg = kset.metadata['symmetry_dataset_pc']
    click.echo('\nPrimitive cell information:')
    click.echo(' ' * 8 + f'Space group number: {pc_spg["number"]}')
    click.echo(' ' * 8 + f'International symbol: {pc_spg["international"]}')
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
