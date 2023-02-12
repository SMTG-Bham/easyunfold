"""
Tests for the CLI system
"""
import os
import pytest
from pathlib import Path

from monty.serialization import loadfn
from click.testing import CliRunner
from easyunfold.cli import easyunfold
from easyunfold.unfold import read_kpoints


def test_generate(si_project_dir):
    """
    Test the generate function
    """
    runner = CliRunner()

    tmp_dir = si_project_dir(None)
    os.chdir(tmp_dir)

    output = runner.invoke(easyunfold, ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test'])
    assert output.exit_code == 0

    output = runner.invoke(
        easyunfold, ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test', '--nk-per-split', '3'])
    assert output.exit_code == 0

    kpts = read_kpoints('KPOINTS_test')[0]
    assert len(kpts) == 11
    kpts = read_kpoints('KPOINTS_test_002')[0]
    assert len(kpts) == 3

    # Test with SCF kpoints

    output = runner.invoke(easyunfold, [
        'generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test', '--scf-kpoints',
        'Si_super_deformed_soc/IBZKPT'
    ])
    assert output.exit_code == 0
    kpts, _, _, weights = read_kpoints('KPOINTS_test')
    assert weights[0] == 12
    assert weights[-1] == 0
    assert len(kpts) == 88 + 11


@pytest.mark.parametrize('tag', ['', '_spin', '_soc'])
def test_unfold(si_project_dir, tag):
    """
    Test the generate function
    """
    runner = CliRunner()

    tmp_dir = si_project_dir(tag)
    os.chdir(tmp_dir)

    output = runner.invoke(easyunfold, ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test.json'])

    # Perform the unfold
    args_calc = ['unfold', '--data-file', 'test.json', 'calculate', f'Si_super_deformed{tag}/WAVECAR']
    if 'soc' in tag:
        args_calc.append('--ncl')
    output = runner.invoke(easyunfold, args_calc)
    assert output.exit_code == 0

    unfoldset = loadfn('test.json')

    assert unfoldset.is_calculated

    # Effective mass
    if tag == '':
        output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'effective-mass'])
        assert 'Hole effective masses' in output.stdout
        assert 'm_e             -1.38922            32  [0.0, 0.0, 0.0] (\Gamma)  [0.5, 0.5, 0.5] (L)' in output.stdout

    # Do the plotting
    output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'plot'])
    assert output.exit_code == 0
    print(output.stdout)
    assert Path('unfold.png').is_file()
