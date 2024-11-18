"""
Tests for the CLI system
"""
import os
from pathlib import Path
import pytest

from monty.serialization import loadfn
from click.testing import CliRunner
from easyunfold.cli import easyunfold
from easyunfold.utils import read_kpoints


def test_generate(si_project_dir):
    """
    Test the generate function
    """
    runner = CliRunner()

    tmp_dir = si_project_dir(None)
    os.chdir(tmp_dir)

    output = runner.invoke(
        easyunfold,
        [
            'generate',
            'Si/POSCAR',
            'Si_super_deformed/POSCAR',
            'KPOINTS_band_low',
            '--out-file',
            'test',
            '-y',
        ],
    )
    assert output.exit_code == 0

    output = runner.invoke(
        easyunfold,
        [
            'generate',
            'Si/POSCAR',
            'Si_super_deformed/POSCAR',
            'KPOINTS_band_low',
            '--out-file',
            'test',
            '--nk-per-split',
            '3',
            '-y',
        ],
    )
    assert output.exit_code == 0

    kpts = read_kpoints('KPOINTS_test')[0]
    kpts_expected = 25
    assert len(kpts) == kpts_expected
    kpts = read_kpoints('KPOINTS_test_002')[0]
    assert len(kpts) == 3

    # Test with SCF kpoints

    output = runner.invoke(
        easyunfold,
        [
            'generate',
            'Si/POSCAR',
            'Si_super_deformed/POSCAR',
            'KPOINTS_band_low',
            '--out-file',
            'test',
            '--scf-kpoints',
            'Si_super_deformed_soc/IBZKPT',
            '-y',
        ],
    )
    assert output.exit_code == 0
    kpts, _, _, weights = read_kpoints('KPOINTS_test')
    assert weights[0] == 1.0
    assert weights[-1] == 0
    # IBZKPT  96 band 28
    assert len(kpts) == 96 + kpts_expected


def test_generate_agsbte2(agsbte2_project_dir):
    """
    Test the generate function with SQS AgSbTe2 – in particular the ability to handle slightly-incommensurate cells
    """
    runner = CliRunner()

    tmp_dir = agsbte2_project_dir(None)
    os.chdir(tmp_dir)

    output = runner.invoke(
        easyunfold,  # ~2.5% lattice mismatch in this case, print warning but continue fine
        ['generate', 'POSCAR_prim', 'SQS_POSCAR', 'KPOINTS_band', '--out-file', 'test'],
    )
    assert output.exit_code == 0

    _kpts = _check_output_info_and_kpoints_agsbte2(
        [
            'Warning: There is a lattice parameter mismatch in the range 2-5% between the primitive (multiplied by the '
            'transformation matrix) and the supercell. This will lead to some quantitative inaccuracies in the '
            'Brillouin Zone spacing (and thus effective masses) of the unfolded band structures.',
            '(Guessed) Transform matrix:\n[[1.0, 0.0, 0.0], [1.0, -3.0, 1.0], [1.0, 1.0, -3.0]]',
        ],
        output,
    )
    # test with also specifying transformation matrix:
    os.remove('KPOINTS_test')
    os.remove('test')
    output = runner.invoke(
        easyunfold,  # ~2.5% lattice mismatch in this case, print warning but continue fine
        [
            'generate',
            'POSCAR_prim',
            'SQS_POSCAR',
            'KPOINTS_band',
            '--out-file',
            'test',
            '-m',
            '1 0 0 1 -3 1 1 1 -3',
        ],
    )
    assert output.exit_code == 0

    assert ('Warning: There is a lattice parameter mismatch in the range 2-5% between the primitive (multiplied by the '
            'transformation matrix) and the supercell. This will lead to some quantitative inaccuracies in the '
            'Brillouin Zone spacing (and thus effective masses) of the unfolded band structures.' in output.output)
    _kpts = _check_output_info_and_kpoints_agsbte2(
        [
            'Proceeding with the assumed transformation matrix.',
            'Transform matrix:\n[[1.0, 0.0, 0.0], [1.0, -3.0, 1.0], [1.0, 1.0, -3.0]]',
        ],
        output,
    )


def _check_agsbte2_info_printing(output):
    assert '194 kpoints specified along the path' in output.output
    for i in [
            'Supercell cell information:',
            'Space group number: 6',
            'International symbol: Pm',
            'Point group: m',
            'Primitive cell information:',
            'Space group number: 225',
            'International symbol: Fm-3m',
            'Point group: m-3m',
    ]:
        assert i in output.output
    assert 'Supercell kpoints written to KPOINTS_test' in output.output
    assert 'Unfolding settings written to test' in output.output

    return read_kpoints('KPOINTS_test')[0]


def _check_output_info_and_kpoints_agsbte2(info_messages: list, output: str):
    for info_message in info_messages:
        assert info_message in output.output
    kpts = _check_agsbte2_info_printing(output)
    assert len(kpts) == 993

    return kpts


def test_generate_from_zero_weighted_prim_kpoints(alvfe_test_dir):
    """
    Test KPOINTS generation when the input KPOINTS are from a meta-GGA
    or hybrid DFT VASP calculation (having weighted and zero-weighted
    kpoints).
    """
    runner = CliRunner()
    tmp_dir = alvfe_test_dir(None)
    os.chdir(tmp_dir)

    output = runner.invoke(
        easyunfold,
        ['generate', 'prim_POSCAR', 'super_POSCAR', 'prim_KPOINTS'],
    )
    kpts_expected = 208
    assert output.exit_code == 0
    assert 'Only zero-weighted kpoints are taken from the input KPOINTS.' in output.output
    assert f'{kpts_expected} kpoints specified along the path' in output.output
    kpts, _, _, weights = read_kpoints('KPOINTS_easyunfold')
    assert len(kpts) == len(weights)
    assert len(kpts) == 324


@pytest.mark.parametrize('tag', ['', '_spin', '_soc'])
def test_unfold(si_project_dir, tag):
    """
    Test the generate and unfolding functions
    """
    runner = CliRunner()

    tmp_dir = si_project_dir(tag)
    os.chdir(tmp_dir)

    output = runner.invoke(
        easyunfold,
        [
            'generate',
            'Si/POSCAR',
            'Si_super_deformed/POSCAR',
            'KPOINTS_band_low',
            '--out-file',
            'test.json',
        ],
    )

    # Status check
    args_calc = ['unfold', '--data-file', 'test.json', 'status']
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0
    assert 'Please run the supercell' in output.stdout

    # Perform the unfold
    args_calc = [
        'unfold',
        '--data-file',
        'test.json',
        'calculate',
        f'Si_super_deformed{tag}/WAVECAR',
    ]
    if 'soc' in tag:
        args_calc.append('--ncl')
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0

    # Status check after unfolding is done
    args_calc = ['unfold', '--data-file', 'test.json', 'status']
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0
    assert 'has been performed' in output.stdout

    unfoldset = loadfn('test.json')

    assert unfoldset.is_calculated

    # Effective mass
    if tag == '':
        output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'effective-mass'])
        assert 'Hole effective masses' in output.stdout
        assert r'0  m_e                0.82              8  [0.5, 0.0, 0.5] (X)' in output.stdout
        # Plot effective mass
        output = runner.invoke(
            easyunfold,
            ['unfold', '--data-file', 'test.json', 'effective-mass', '--plot'],
        )
        assert Path('unfold-effective-mass.png').is_file()
        # Plot effective mass fits
        output = runner.invoke(
            easyunfold,
            ['unfold', '--data-file', 'test.json', 'effective-mass', '--plot-fit'],
        )
        assert Path('unfold-effective-mass.png').is_file()

    # Do the plotting
    output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'plot'])
    assert output.exit_code == 0
    print(output.stdout)
    assert Path('unfold.png').is_file()

    # matplotlib customisation check
    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'test.json',
            'plot',
            '--mpl-style-file',
            'my.mplstyle',
        ],
    )
    assert output.exit_code == 0
    print(output.stdout)
    assert 'Using custom plotting style' in output.stdout
    assert Path('unfold.png').is_file()


def test_plot_projection(mgo_project_dir):
    """Test plot projection"""
    os.chdir(mgo_project_dir)
    runner = CliRunner()
    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'mgo.json',
            'plot-projections',
            '--atoms-idx',
            '1,2|3-4',
            '--procar',
            'PROCAR',
        ],
    )
    _plot_projection_check(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'mgo.json',
            'plot-projections',
            '--atoms-idx',
            '1,2|3-4',
            '--procar',
            'PROCAR',
            '--combined',
        ],
    )
    _plot_projection_check(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'mgo.json',
            'plot-projections',
            '--atoms-idx',
            '1,2|3-4',
            '--procar',
            'PROCAR',
            '--combined',
            '--orbitals',
            'px,py|pz',
        ],
    )
    _plot_projection_check(output)

    # test --atoms option with --poscar specification
    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'mgo.json',
            'plot-projections',
            '--atoms',
            'Mg,O',
            '--poscar',
            'POSCAR.mgo',
        ],
    )
    _plot_projection_check(output)

    # test parsing PROCAR from LORBIT = 14 calculation
    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'mgo.json',
            'plot-projections',
            '--atoms',
            'Mg,O',
            '--poscar',
            'POSCAR.mgo',
            '--procar',
            'PROCAR_LORBIT_14.mgo',
        ],
    )
    _plot_projection_check(output)

    # test options
    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            '--data-file',
            'mgo.json',
            'plot-projections',
            '--atoms',
            'Mg,O',
            '--poscar',
            'POSCAR.mgo',
            '--eref',
            '2',
            '--no-symm-average',
            '--cmap',
            'PuBu',
            '--orbitals',
            's,p',
            '--colours',
            'r,g',
            '--orbitals',
            's,p',
            '--colours',
            'r,g',
            '--colourspace',
            'luvlch',
            '--intensity',
            '0.5',
            '--emin',
            '-2',
            '--emax',
            '5',
            '--dpi',
            '500',
            '--vscale',
            '2.0',
            '--cmap',
            'PuBu',
            '--npoints',
            '200',
            '--sigma',
            '0.15',
            '--title',
            'Test',
            '--no-combined',
            '--height',
            '2',
            '--width',
            '3',
            '-o',
            'test.png',
        ],
    )
    assert output.exit_code == 0
    assert Path('test.png').is_file()  # different file name this time
    Path('test.png').unlink()


def _plot_projection_check(output):
    assert output.exit_code == 0
    assert Path('unfold.png').is_file()
    Path('unfold.png').unlink()


def test_help(nabis2_project_dir):
    """Test help messages"""
    os.chdir(nabis2_project_dir)
    runner = CliRunner()
    output = runner.invoke(easyunfold, ['unfold', 'plot-projections', '-h'])
    assert output.exit_code == 0
    assert not Path('unfold.png').is_file()

    output = runner.invoke(easyunfold, ['unfold', 'plot', '-h'])
    assert output.exit_code == 0
    assert not Path('unfold.png').is_file()

    output = runner.invoke(easyunfold, ['unfold', 'plot-projections', '--help'])
    assert output.exit_code == 0
    assert not Path('unfold.png').is_file()

    output = runner.invoke(easyunfold, ['unfold', 'plot', '--help'])
    assert output.exit_code == 0
    assert not Path('unfold.png').is_file()


def test_dos_atom_orbital_plots(nabis2_project_dir):
    """Test various dos/atom/orbital etc plot options with NaBiS2"""
    os.chdir(nabis2_project_dir)
    runner = CliRunner()
    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi,S',
            '--orbitals',
            's|px,py,pz|p',
            '--vscale',
            '0.5',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
            '--zero-line',
            '--dos-label',
            'DOS',
            '--gaussian',
            '0.1',
            '--no-total',
            '--scale',
            '2',
            '--dos-elements',
            'Bi.s.p',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi,S',
            '--orbitals',
            's|px,py,pz|p',
            '--vscale',
            '0.5',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
            '--zero-line',
            '--dos-label',
            'DOS',
            '--gaussian',
            '0.1',
            '--no-total',
            '--scale',
            '2',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [  # same but with intensity instead of vscale:
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi,S',
            '--orbitals',
            's|px,py,pz|p',
            '--intensity',
            '2',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
            '--zero-line',
            '--dos-label',
            'DOS',
            '--gaussian',
            '0.1',
            '--no-total',
            '--scale',
            '2',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi,S',
            '--vscale',
            '0.5',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
            '--zero-line',
            '--dos-label',
            'DOS',
            '--gaussian',
            '0.1',
            '--no-total',
            '--scale',
            '2',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi',
            '--vscale',
            '0.5',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
            '--zero-line',
            '--dos-label',
            'DOS',
            '--gaussian',
            '0.1',
            '--no-total',
            '--scale',
            '2',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        ['unfold', 'plot-projections', '--atoms', 'Na,Bi,S', '--dos', 'vasprun.xml.gz'],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi,S',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot',
            '--atoms-idx',
            '1-20|21-40',
            '--orbitals',
            's|p',
            '--dos',
            'vasprun.xml.gz',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot',
            '--atoms',
            'Na,Bi',
            '--orbitals',
            's|p',
            '--dos',
            'vasprun.xml.gz',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi',
            '--combined',
            '--orbitals',
            's',
            '--dos',
            'vasprun.xml.gz',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi',
            '--combined',
            '--orbitals',
            's',
            '--dos',
            'vasprun.xml.gz',
            '--dos-elements',
            'Bi.s',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms-idx',
            '1-20,21,22,33',
            '--combined',
            '--orbitals',
            's',
            '--dos',
            'vasprun.xml.gz',
            '--dos-elements',
            'Bi.s',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot',
            '--atoms-idx',
            '1-20,21,22,33',
            '--orbitals',
            's',
            '--dos',
            'vasprun.xml.gz',
            '--dos-elements',
            'Bi.s',
        ],
    )
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(easyunfold, ['unfold', 'plot', '--dos', 'vasprun.xml.gz'])
    _check_dos_atom_orbital_plots(output)

    output = runner.invoke(
        easyunfold,
        [
            'unfold',
            'plot-projections',
            '--atoms',
            'Na,Bi',
            '--orbitals',
            's',
            '--combined',
            '--dos',
            'vasprun.xml.gz',
            '--dos-elements',
            'Bi.s',
        ],
    )
    _check_dos_atom_orbital_plots(output)


def _check_dos_atom_orbital_plots(output):
    assert output.exit_code == 0
    assert Path('unfold.png').is_file()
    Path('unfold.png').unlink()
