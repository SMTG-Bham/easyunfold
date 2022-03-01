"""
Tests for the CLI system
"""
import os
import urllib
import pytest
from pathlib import Path
import shutil

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

    kpts, _, _ = read_kpoints('KPOINTS_test')
    assert len(kpts) == 11


@pytest.mark.parametrize('tag', ['', '_spin', '_soc'])
def test_generate(si_project_dir, tag):
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

    # Do the plotting
    output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'plot'])
    assert output.exit_code == 0
    print(output.stdout)
    assert Path('unfold.png').is_file()
