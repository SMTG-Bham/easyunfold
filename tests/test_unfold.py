"""
Test unfolding routines
"""
import shutil
import urllib.request
import numpy as np
from ase.io import read
import pytest
import easyunfold.unfold as unfold


@pytest.fixture
def si_atoms(datapath):
    return read(datapath('Si/POSCAR'))


@pytest.fixture
def si222_atoms(si_atoms):
    return si_atoms.repeat((2, 2, 2))


@pytest.fixture
def si_project_dir(datapath, tmp_path):
    """Create a temporary directory containing the Si unfold project data"""
    shutil.copytree(datapath('Si-project'), tmp_path / 'Si-project')
    if not list(tmp_path.glob('Si-project/*/WAVECAR')):
        # Download the dataset
        data_depo = [
            ('https://www.dropbox.com/s/0d6cc8rsee2j7to/WAVCAR?dl=1', 'Si_super_deformed/WAVECAR'),
            ('https://www.dropbox.com/s/22u33579kf3zq4x/WAVECAR?dl=1', 'Si_super_deformed_spin/WAVECAR'),
        ]
        for url, relpath in data_depo:
            print(f'Downloading {relpath}')
            urllib.request.urlretrieve(url, tmp_path / 'Si-project' / relpath)
    return tmp_path / 'Si-project'


@pytest.fixture
def kpath_and_labels():
    path = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 0.25, 0.75], [0.5, 0.0, 0.5]]
    labels = ['G', 'L', 'W', 'X']
    return path, labels


@pytest.fixture
def explicit_kpoints(kpath_and_labels):
    """Test generating kpoints"""
    klist, klabel = kpath_and_labels
    kpts = unfold.make_kpath(klist, 20)
    assert len(kpts) == 61
    return kpts


@pytest.fixture
def explicit_kpoints_minimal(kpath_and_labels):
    """Test generating kpoints"""
    klist, klabel = kpath_and_labels
    kpts = unfold.make_kpath(klist, 2)
    assert len(kpts) == 7
    return kpts


@pytest.fixture
def silicon_unfold(explicit_kpoints_minimal, si_atoms, si222_atoms):
    """
    Return an object for unfolding silicon
    """
    return unfold.UnfoldKSet.from_atoms(np.diag([2, 2, 2]), explicit_kpoints_minimal, si_atoms, si222_atoms)


def test_read_kpoints_line(datapath):

    kpoints, comment, labels = unfold.read_kpoints(datapath('KPOINTS_LINE'))
    assert len(kpoints) == 91
    assert len(labels) == 4
    assert labels[0][0] == 0
    assert labels[-1][0] == 90
    assert comment == 'Monkhorst-Pack'


def test_unfold_expansion(si_atoms, si222_atoms, explicit_kpoints):
    """Test genearting extended kpoints set"""
    # Get symmetry operations
    rots_pc = unfold.get_symmetry_dataset(si_atoms)['rotations']
    # Break the symmetry
    si222_atoms[0].position += np.array([0.1, 0.1, 0.1])
    rots_sc = unfold.get_symmetry_dataset(si222_atoms)['rotations']

    foldset = unfold.UnfoldKSet(np.diag((2, 2, 2)), explicit_kpoints, si222_atoms.cell, rots_pc, rots_sc)
    assert rots_pc.shape[0] == 48
    assert rots_sc.shape[0] == 6
    assert foldset.nkpts_orig == len(explicit_kpoints)
    assert foldset.nkpts_expand == 138


def test_symmetry_expand(si_atoms, si222_atoms):
    """Test the expansion of symmetrically related points"""
    rots_pc = unfold.get_symmetry_dataset(si_atoms)['rotations']
    rots_sc = unfold.get_symmetry_dataset(si222_atoms)['rotations']
    kpts, weights = unfold.expand_K_by_symmetry([0.1, 0.1, 0.1], rots_pc, rots_sc, time_reversal=True)
    assert len(kpts) == 1
    assert len(weights) == 1

    si222_atoms[0].position += np.array([0.1, 0.1, 0.1])
    rots_pc = unfold.get_symmetry_dataset(si_atoms)['rotations']
    rots_sc = unfold.get_symmetry_dataset(si222_atoms)['rotations']
    kpts, weights = unfold.expand_K_by_symmetry([0.1, 0.1, 0.1], rots_pc, rots_sc, time_reversal=True)
    assert len(kpts) == 2
    assert len(weights) == 2

    kpts, weights = unfold.expand_K_by_symmetry([0.1, 0.1, 0.1], rots_pc, rots_sc, time_reversal=False)
    assert len(kpts) == 4
    assert len(weights) == 4


def test_serialization(silicon_unfold, tmp_path):
    """
    Test serializing and loading from the archive file
    """
    from monty.serialization import loadfn
    assert silicon_unfold.as_dict()
    new_obj = silicon_unfold.from_dict(silicon_unfold.as_dict())

    # Test writing out
    (tmp_path / 'out.json').write_text(silicon_unfold.to_json())

    new_obj = loadfn(tmp_path / 'out.json')
    np.testing.assert_allclose(new_obj.M, silicon_unfold.M)
    assert 'kpoints' in new_obj.expansion_results


@pytest.mark.parametrize('folder_name,nspin', [('Si_super_deformed', 1), ('Si_super_deformed_spin', 2)])
def test_unfold(si_project_dir, folder_name, nspin):
    """
    Test unfolding on the real data
    """

    atoms_primitive = read(si_project_dir / 'Si/POSCAR')
    atoms_supercell = read(si_project_dir / f'{folder_name}/POSCAR')
    kpoints, _, labels = unfold.read_kpoints(si_project_dir / 'KPOINTS_band_low')

    unfolder: unfold.UnfoldKSet = unfold.UnfoldKSet.from_atoms(np.diag([2, 2, 2]), kpoints, atoms_primitive, atoms_supercell)
    unfolder.kpoint_labels = labels
    unfolder.write_sc_kpoints(si_project_dir / 'KPOINTS_sc')

    # Test kpoints generation
    kpoints_sc, _, _ = unfold.read_kpoints(si_project_dir / 'KPOINTS_sc')
    kpoints_sc_ref, _, _ = unfold.read_kpoints(si_project_dir / f'{folder_name}/KPOINTS_easyunfold')
    np.testing.assert_allclose(kpoints_sc, kpoints_sc_ref)

    # Test unfold
    sws = unfolder.get_spectral_weights(si_project_dir / f'{folder_name}/WAVECAR')
    assert sws.shape[0] == nspin
    assert sws.shape[1] == len(kpoints)
    assert sws.shape[2] == 81
    assert sws.shape[3] == 2

    assert unfolder.is_calculated

    # Spectral weights
    e0, spectral_function = unfolder.get_spectral_function(npoints=500)
    assert len(e0) == 500
    assert spectral_function.shape == (nspin, 500, len(kpoints))

    assert unfolder.to_json()
