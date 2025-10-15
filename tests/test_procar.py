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

    assert nb == 20
    assert nk == 47

    assert procar.nion == 2
    assert procar.eigenvalues.shape == (1, nk, nb)
    assert procar.kvecs.shape == (1, nk, 3)
    assert procar.kweights.shape == (1, nk)
    assert np.all(procar.kvecs[0][0] == 0.0)
    assert procar.occs.shape == (1, nk, nb)
    assert procar.proj_data[0, 0, 0, 0, 0] == 0.041
    assert procar.proj_data[0, 0, 0, 1, 0] == 0.678

    assert (procar.get_projection([0], 'all').sum() != 618.2850603903657)  # doesn't match before normalisation

    # After normalisation
    procar.normalise_projs(procar.proj_data)
    for nspin in range(procar.proj_data.shape[0]):
        for nband in range(procar.proj_data.shape[1]):
            for nkpt in range(procar.proj_data.shape[2]):
                # either zero projection or normalised to 1:
                assert any(np.isclose(procar.proj_data[nspin, nband, nkpt, :, :].sum(), x) for x in [0, 1])

    assert procar.get_projection([0], 'all').sum() == 618.2850603903657
    assert procar.get_projection([0], ['s', 'px']).sum() == 124.10684519359549
    # assert individual s and px projections are same as sum above with both:
    assert (procar.get_projection([0], ['s']).sum() + procar.get_projection([0], 'px').sum() == 124.10684519359546)
    assert procar.proj_names == [
        's',
        'py',
        'pz',
        'px',
        'dxy',
        'dyz',
        'dz2',
        'dxz',
        'x2-y2',
    ]
    assert not procar._is_soc

    # check projections behave as expected:
    assert np.isclose(
        procar.get_projection(list(range(procar.nion)), 'all').sum(),
        np.sum(procar.proj_data),
    )

    # get number of procar.proj_data entries where [x, y, z, :, :].sum() == 0:  (i.e. zero total projection on all
    # atoms, happens for some bands (3 in this case) at high-symmetry kpoints due to orthogonal projections on
    # spherical harmonics)
    num_zero = np.sum([
        1 for x in range(procar.proj_data.shape[0]) for y in range(procar.proj_data.shape[1]) for z in range(procar.proj_data.shape[2])
        if np.isclose(procar.proj_data[x, y, z, :, :].sum(), 0)
    ])
    assert np.isclose(
        np.sum(procar.proj_data) + num_zero,
        procar.nbands * procar.nkpts * procar.nspins,
    )

    # test weight_by_k:
    assert np.isclose(
        procar.get_projection(list(range(procar.nion)), 'all', weight_by_k=True).sum() *
        (48 / 47)  # because Gamma duplicated in this calc, so reduced to 47 of orig 48 kpoints. Either way, weights
        # not meaningful in this case as it's a non-SCF BS calc, just testing weighting runs as expected
        + num_zero / procar.nkpts,
        procar.nbands * procar.nspins,
    )


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
    assert np.all(procar.kvecs[0][0] == 0.0)
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
    assert np.all(procar.kvecs[0][0] == 0.0)
    assert procar.occs.shape == (2, nk, nb)
