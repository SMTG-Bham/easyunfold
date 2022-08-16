"""
Test the effective mass functionality
"""
from pathlib import Path
import pytest
from monty.serialization import loadfn
import easyunfold.effective_mass as em


def test_fit_effective_mass():
    """Test effective mass fitting"""

    dists = [0, 1, 2]
    eng = [0, 2, 4]
    assert em.fit_effective_mass(dists, eng)


@pytest.fixture(scope='module')
def effective_mass_obj():
    """Return an unfolding object"""
    obj = loadfn(Path(__file__).parent / 'test_data/mgo.json')
    efm = em.EffectiveMass(obj)
    return efm


def test_band_maxima(effective_mass_obj):
    """Test detection of band extreme"""
    output = effective_mass_obj.get_band_extrema()
    assert output == ([0, 47], [0, 0], [16, 16])
    output = effective_mass_obj.get_band_extrema(mode='vbm')
    assert output == ([0, 47], [0, 0], [15, 15])


def test_effective_mass(effective_mass_obj):
    """Test obtaining effective masses"""

    output = effective_mass_obj.get_effective_masses()

    kdist = effective_mass_obj._get_kpoint_distances()
    assert kdist[0] == 0
    assert kdist[1] == pytest.approx(0.09833682)
    assert len(kdist) == len(effective_mass_obj.kpoints)

    fdata = effective_mass_obj._get_fitting_data(0, 16, 1, 0, 3)
    assert lenfdata[0] == kdist[:3]
    assert len(fdata[1]) == 3

    assert len(output['electrons']) == 2
    assert len(output['holes']) == 2

    assert output['electrons'][0]['kpoint_label_from'] == '\\Gamma'
    assert output['electrons'][0]['kpoint_label_to'] == 'X'
    assert output['electrons'][0]['effective_mass'] == 0.3871366744253859

    assert output['holes'][0]['kpoint_label_from'] == '\\Gamma'
    assert output['holes'][0]['kpoint_label_to'] == 'X'
    assert output['holes'][0]['effective_mass'] == -2.8972969334826235
