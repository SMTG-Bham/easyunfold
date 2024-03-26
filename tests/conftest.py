"""
Pytest fixtures
"""

import urllib.request
from pathlib import Path
import shutil

import pytest


@pytest.fixture
def datapath():
    """Return the path relative to the test_data folder"""

    def _inner(rel):

        return Path(__file__).parent / 'test_data' / rel

    return _inner


@pytest.fixture
def si_project_dir(datapath, tmp_path):
    """Create a temporary directory containing the Si unfold project data"""

    def _inner(tag=''):
        shutil.copytree(datapath('Si-project'), tmp_path / 'Si-project')
        # if not (tmp_path / f'Si-project/Si_super_deformed{tag}/WAVECAR').is_file() and tag is not None:
        #     # Download the dataset
        #     relpath = f'Si_super_deformed{tag}/WAVECAR'
        #     url = DATA_REPO[relpath]
        #     print(f'Downloading {relpath}')
        #     urllib.request.urlretrieve(url, tmp_path / 'Si-project' / relpath)
        #     # Copy the downloaded WAVECAR to the repository for future use
        #     shutil.copy(tmp_path / 'Si-project' / relpath, datapath('Si-project') / f'Si_super_deformed{tag}/WAVECAR')
        return tmp_path / 'Si-project'

    return _inner


@pytest.fixture
def agsbte2_project_dir(datapath, tmp_path):
    """Create a temporary directory containing AgSbTe2 inputs for testing"""

    def _inner(tag=''):
        shutil.copytree(datapath('AgSbTe2'), tmp_path / 'AgSbTe2')
        return tmp_path / 'AgSbTe2'

    return _inner


@pytest.fixture
def mgo_project_dir(datapath, tmp_path):
    shutil.copy2(datapath('mgo.json'), tmp_path / 'mgo.json')
    shutil.copy2(datapath('PROCAR.mgo.gz'), tmp_path / 'PROCAR.gz')
    shutil.copy2(datapath('POSCAR.mgo'), tmp_path / 'POSCAR.mgo')
    shutil.copy2(datapath('PROCAR_LORBIT_14.mgo.gz'), tmp_path / 'PROCAR_LORBIT_14.mgo.gz')
    return tmp_path


@pytest.fixture
def nabis2_project_dir(tmp_path):
    nabis2_dir = Path(__file__).parent / '..' / 'examples' / 'NaBiS2'
    for i in ['POSCAR', 'KPOINTS', 'easyunfold.json', 'PROCAR.gz', 'vasprun.xml.gz']:
        shutil.copy2(nabis2_dir / i, tmp_path / i)
    return tmp_path
