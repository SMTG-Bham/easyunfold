"""
Pytest fixtures
"""

import urllib.request
from pathlib import Path
import shutil

import pytest

_data_depo = [
    ('https://www.dropbox.com/s/0d6cc8rsee2j7to/WAVCAR?dl=1', 'Si_super_deformed/WAVECAR'),
    ('https://www.dropbox.com/s/22u33579kf3zq4x/WAVECAR?dl=1', 'Si_super_deformed_spin/WAVECAR'),
    ('https://www.dropbox.com/s/4y271w81bkf0gvx/WAVECAR?dl=1', 'Si_super_deformed_soc/WAVECAR'),
]
DATA_REPO = {key: value for value, key in _data_depo}


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
        if not (tmp_path / f'Si-project/Si_super_deformed{tag}/WAVECAR').is_file() and tag is not None:
            # Download the dataset
            relpath = f'Si_super_deformed{tag}/WAVECAR'
            url = DATA_REPO[relpath]
            print(f'Downloading {relpath}')
            urllib.request.urlretrieve(url, tmp_path / 'Si-project' / relpath)
        return tmp_path / 'Si-project'

    return _inner
