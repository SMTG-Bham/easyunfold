"""
Tests for the CLI system
"""
import os
import pytest
from pathlib import Path

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

    output = runner.invoke(easyunfold,
                           ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test', '-y'])
    assert output.exit_code == 0

    output = runner.invoke(
        easyunfold,
        ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test', '--nk-per-split', '3', '-y'])
    assert output.exit_code == 0

    kpts = read_kpoints('KPOINTS_test')[0]
    kpts_expected = 25
    assert len(kpts) == kpts_expected
    kpts = read_kpoints('KPOINTS_test_002')[0]
    assert len(kpts) == 3

    # Test with SCF kpoints

    output = runner.invoke(easyunfold, [
        'generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test', '--scf-kpoints',
        'Si_super_deformed_soc/IBZKPT', '-y'
    ])
    assert output.exit_code == 0
    kpts, _, _, weights = read_kpoints('KPOINTS_test')
    assert weights[0] == 1.0
    assert weights[-1] == 0
    # IBZKPT  96 band 28
    assert len(kpts) == 96 + kpts_expected


@pytest.mark.parametrize('tag', ['', '_spin', '_soc'])
def test_unfold(si_project_dir, tag):
    """
    Test the generate function
    """
    runner = CliRunner()

    tmp_dir = si_project_dir(tag)
    os.chdir(tmp_dir)

    output = runner.invoke(easyunfold, ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test.json'])

    # Status check
    args_calc = ['unfold', '--data-file', 'test.json', 'status']
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0
    assert 'Please run the supercell' in output.stdout

    # Perform the unfold
    args_calc = ['unfold', '--data-file', 'test.json', 'calculate', f'Si_super_deformed{tag}/WAVECAR']
    if 'soc' in tag:
        args_calc.append('--ncl')
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0

    # Status check after unfolding is done
    args_calc = ['unfold', '--data-file', 'test.json', 'status']
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0
    assert 'had been performed' in output.stdout

    unfoldset = loadfn('test.json')

    assert unfoldset.is_calculated

    # Effective mass
    if tag == '':
        output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'effective-mass'])
        assert 'Hole effective masses' in output.stdout
        assert r'     0  m_e            -0.938459             8  [0.5, 0.0, 0.5] (X)  [0.5, 0.25, 0.75] (W)' in output.stdout
        # Plot effective mass
        output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'effective-mass', '--plot'])
        assert Path('unfold-effective-mass.png').is_file()
        # Plot effective mass fits
        output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'effective-mass', '--plot-fit'])
        assert Path('unfold-effective-mass.png').is_file()

    # Do the plotting
    output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'plot'])
    assert output.exit_code == 0
    print(output.stdout)
    assert Path('unfold.png').is_file()


def test_plot_projection(mgo_project_dir):
    """Test plot projection"""
    os.chdir(mgo_project_dir)
    runner = CliRunner()
    output = runner.invoke(easyunfold,
                           ['unfold', '--data-file', 'mgo.json', 'plot-projections', '--atoms-idx', '1,2|3-4', '--procar', 'PROCAR'])
    assert output.exit_code == 0
    assert Path('unfold.png').is_file()
    Path('unfold.png').unlink()

    output = runner.invoke(
        easyunfold, ['unfold', '--data-file', 'mgo.json', 'plot-projections', '--atoms-idx', '1,2|3-4', '--procar', 'PROCAR', '--combined'])
    assert output.exit_code == 0
    assert Path('unfold.png').is_file()
    Path('unfold.png').unlink()

    output = runner.invoke(easyunfold, [
        'unfold', '--data-file', 'mgo.json', 'plot-projections', '--atoms-idx', '1,2|3-4', '--procar', 'PROCAR', '--combined', '--orbitals',
        'px,py|pz'
    ])
    assert output.exit_code == 0
    assert Path('unfold.png').is_file()
    Path('unfold.png').unlink()
