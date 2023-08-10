"""
Tests from reading PROCAR
"""

from pathlib import Path
import numpy as np

from easyunfold.procar import Procar

datapath = Path(__file__).parent / 'test_data'


def test_procar():
    """Test reading PROCAR file"""

    procar = Procar(datapath / 'PROCAR')

    assert procar.nion == 2
    assert procar.eigenvalues.shape == (1, 47, 20)
    assert procar.kvecs.shape == (1, 47, 3)
    assert procar.kweights.shape == (1, 47)
    assert np.all(procar.kvecs[0][0] == 0.)
    assert procar.occs.shape == (1, 47, 20)
    assert procar.get_projection([0], 'all').sum() == 618.2850603903657
    assert procar.get_projection([0], ['s', 'px']).sum() == 124.10684519359549
    assert procar.get_projection([0], ['s']).sum() + procar.get_projection([0], 'px').sum() == 124.10684519359546
    assert procar.proj_names == ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
