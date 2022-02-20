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


@pytest.fixture
def si_project_dir(datapath, tmp_path):
    """Create a temporary directory containing the Si unfold project data"""

    def _run(download=False):
        shutil.copytree(datapath('Si-project'), tmp_path / 'Si-project')
        if not list(tmp_path.glob('Si-project/*/WAVECAR')) and download:
            # Download the dataset
            data_depo = [
                ('https://www.dropbox.com/s/0d6cc8rsee2j7to/WAVCAR?dl=1', 'Si_super_deformed/WAVECAR'),
                ('https://www.dropbox.com/s/22u33579kf3zq4x/WAVECAR?dl=1', 'Si_super_deformed_spin/WAVECAR'),
            ]
            for url, relpath in data_depo:
                print(f'Downloading {relpath}')
                urllib.request.urlretrieve(url, tmp_path / 'Si-project' / relpath)
        return tmp_path / 'Si-project'

    return _run


def test_generate(si_project_dir):
    """
    Test the generate function
    """
    runner = CliRunner()

    tmp_dir = si_project_dir()
    os.chdir(tmp_dir)

    output = runner.invoke(easyunfold, ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test'])
    assert output.return_value is None

    kpts, _, _ = read_kpoints('KPOINTS_test')
    assert len(kpts) == 11


def test_generate(si_project_dir):
    """
    Test the generate function
    """
    runner = CliRunner()

    tmp_dir = si_project_dir(True)
    os.chdir(tmp_dir)

    output = runner.invoke(easyunfold, ['generate', 'Si/POSCAR', 'Si_super_deformed/POSCAR', 'KPOINTS_band_low', '--out-file', 'test.json'])

    # Perform the unfold
    output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'calculate', 'Si_super_deformed/WAVECAR'])
    assert output.return_value is None

    unfoldset = loadfn('test.json')

    assert unfoldset.is_calculated

    # Do the plotting
    output = runner.invoke(easyunfold, ['unfold', '--data-file', 'test.json', 'plot'])
    assert output.return_value is None
    print(output.stdout)
    assert Path('unfold.png').is_file()
