"""
Tests from reading PROCAR
"""

from pathlib import Path
import numpy as np

from easyunfold.procar import Procar

datapath = Path(__file__).parent / 'test_data'


def test_procar():
    """Test reading PROCAR file"""

    procar = Procar(datapath / 'PROCAR', normalise=False)

    nb = procar.nbands
    nk = procar.nkpts

    assert procar.nion == 2
    assert procar.eigenvalues.shape == (1, nk, nb)
    assert procar.kvecs.shape == (1, nk, 3)
    assert procar.kweights.shape == (1, nk)
    assert np.all(procar.kvecs[0][0] == 0.)
    assert procar.occs.shape == (1, nk, nb)
    assert procar.proj_data[0, 0, 0, 0, 0] == 0.041
    assert procar.proj_data[0, 0, 0, 1, 0] == 0.678

    # After normalisation
    procar.normalise_projs(procar.proj_data)
    assert procar.get_projection([0], 'all').sum() == 629.7023061524948
    assert procar.get_projection([0], ['s', 'px']).sum() == 126.61032305675073
    assert procar.get_projection([0], ['s']).sum() + procar.get_projection([0], 'px').sum() == 63.92486489016277 + 62.685458166587956
    assert procar.proj_names == ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
    assert not procar._is_soc


def test_procar_soc():
    """Test reading PROCAR file"""

    procar = Procar(datapath / 'Si-project/Si_super_deformed_soc/PROCAR')

    assert procar.nion == 4
    assert procar.nkpts == 25
    assert procar.nspins == 1
    assert procar._is_soc
    nb = procar.nbands
    nk = procar.nkpts
    assert procar.eigenvalues.shape == (1, nk, nb)
    assert procar.kvecs.shape == (1, nk, 3)
    assert procar.kweights.shape == (1, nk)
    assert np.all(procar.kvecs[0][0] == 0.)
    assert procar.occs.shape == (1, nk, nb)


def test_procar_spin():
    """Test reading PROCAR file"""

    procar = Procar(datapath / 'Si-project/Si_super_deformed_spin/PROCAR')

    assert procar.nion == 4
    assert procar.nkpts == 25
    assert procar.nspins == 2
    assert not procar._is_soc
    nb = procar.nbands
    nk = procar.nkpts
    assert procar.eigenvalues.shape == (2, nk, nb)
    assert procar.kvecs.shape == (2, nk, 3)
    assert procar.kweights.shape == (2, nk)
    assert np.all(procar.kvecs[0][0] == 0.)
    assert procar.occs.shape == (2, nk, nb)
