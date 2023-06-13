"""
Tests for utility module
"""

from easyunfold import utils
import numpy as np


def test_find_unique():
    """Test for the find unique function"""
    data = np.random.rand(100, 2)
    data = np.round(data, 3)

    unique, idx, idx_inv = utils.find_unique(data, lambda x, y: np.all(x == y))

    assert np.all(np.in1d(np.unique(data, axis=0), unique))
    assert np.all(np.in1d(unique, np.unique(data, axis=0)))
    assert len(np.unique(data, axis=0)) == len(unique)
    assert np.all(data == unique[idx_inv])
    assert np.all(unique == data[idx])


def test_wrap_kpoints():
    """Test wrapping kpoints"""
    kpoints = np.random.rand(100, 3)
    wrapped = utils.wrap_kpoints(kpoints)
    assert wrapped.max() < 0.5
    assert wrapped.min() >= -0.5
    # Wrap both in range 0, 1
    diff = ((wrapped - np.floor(wrapped)) - (kpoints - np.floor(kpoints)))
    # Make sure they all agree with each other
    assert abs(diff).max() < 1e-10


def test_reduce_kpoints():
    """Test kpoints reduction via wrapping + time-reversal"""
    kpoints = np.array([[0, 0, 0], [0.1, 0.1, 0.1], [-0.1, -0.1, -0.1], [1.1, 1.1, 1.1]])
    reduced, unique_id, inv_mapping = utils.reduce_kpoints(kpoints)
    assert len(reduced) == 2
    assert np.all(unique_id == np.array([0, 1]))
    assert np.all(inv_mapping == [0, 1, 1, 1])

    # Note that np.unique will sort the values so the order of the returned objects are different
    # than that of utils.find_unique
    reduced, unique_id, inv_mapping = utils.reduce_kpoints(kpoints, time_reversal=False)
    assert len(reduced) == 3
    assert 0 in unique_id
    assert 1 in unique_id
    assert 2 in unique_id
