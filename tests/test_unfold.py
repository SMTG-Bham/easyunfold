"""
Test unfolding routines
"""
import numpy as np
from ase.io import read
import pytest
import shutil
import easyunfold.unfold as unfold
from easyunfold.utils import read_kpoints
from matplotlib.colors import hex2color


@pytest.fixture
def si_atoms(datapath):
    return read(datapath('Si/POSCAR'))


@pytest.fixture
def si222_atoms(si_atoms):
    return si_atoms.repeat((2, 2, 2))


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
    """Test reading kpoints in the line mode"""
    kpoints, comment, labels, _ = read_kpoints(datapath('KPOINTS_LINE'))
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
    assert foldset.nkpts_expand == 119


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


@pytest.mark.parametrize('tag,nspin,ncl, nbands_expected', [('', 1, False, 14), ('_spin', 2, False, 14), ('_soc', 1, True, 20)])
def test_unfold(si_project_dir, tag, nspin, ncl, nbands_expected, datapath):
    """
    Test unfolding on the real data
    """
    si_project_dir = si_project_dir(tag)
    folder_name = f'Si_super_deformed{tag}'

    atoms_primitive = read(si_project_dir / 'Si/POSCAR')
    atoms_supercell = read(si_project_dir / f'{folder_name}/POSCAR')
    kpoints, _, labels, _ = read_kpoints(si_project_dir / 'KPOINTS_band_low')

    unfolder: unfold.UnfoldKSet = unfold.UnfoldKSet.from_atoms(np.diag([2, 1, 1]), kpoints, atoms_primitive, atoms_supercell)
    unfolder.kpoint_labels = labels
    unfolder.write_sc_kpoints(si_project_dir / 'KPOINTS_sc')
    # Test split
    unfolder.write_sc_kpoints(si_project_dir / 'KPOINTS_sc', nk_per_split=3)
    assert (si_project_dir / 'KPOINTS_sc_001').is_file()
    assert (si_project_dir / 'KPOINTS_sc_002').is_file()
    ktmp1 = read_kpoints(si_project_dir / 'KPOINTS_sc_001')[0]
    assert len(ktmp1) == 3

    # Split with SCF kpoints
    unfolder.write_sc_kpoints(si_project_dir / 'KPOINTS_sc',
                              nk_per_split=3,
                              scf_kpoints_and_weights=([[0., 0., 0.], [0.1, 0.1, 0.1]], [1, 2]))
    ktmp1 = read_kpoints(si_project_dir / 'KPOINTS_sc_001')[0]
    assert len(ktmp1) == 5
    np.testing.assert_allclose(ktmp1[1], [0.1, 0.1, 0.1])

    # Test kpoints generation
    kpoints_sc = read_kpoints(si_project_dir / 'KPOINTS_sc')[0]
    # to update:
    # shutil.copyfile(si_project_dir / 'KPOINTS_sc',
    #                 datapath('Si-project') / f'{folder_name}/KPOINTS_easyunfold')
    kpoints_sc_ref = read_kpoints(si_project_dir / f'{folder_name}/KPOINTS_easyunfold')[0]
    for kpt in kpoints_sc_ref:
        found = False
        for ref in kpoints_sc:
            if np.allclose(kpt, ref, atol=1e-9):
                found = True
                break
        assert found, f'Kpoint {kpt} not found in the unfolded kpoints'

    # Test unfold
    sws = unfolder.get_spectral_weights(si_project_dir / f'{folder_name}/WAVECAR', ncl=ncl)
    assert sws[0].shape[0] == nspin
    assert len(sws) == len(kpoints)
    assert sws[0].shape[2] == nbands_expected
    assert sws[0].shape[3] == 2

    assert unfolder.is_calculated

    # Spectral weights
    e0, spectral_function = unfolder.get_spectral_function(npoints=500, ncl=ncl)
    assert len(e0) == 500
    assert spectral_function.shape == (nspin, 500, len(kpoints))

    assert unfolder.to_json()


@pytest.mark.parametrize('tag,nspin,ncl, nbands_expected', [('_castep', 1, False, 14)])
def test_unfold_castep(si_project_dir, tag, nspin, ncl, nbands_expected, datapath):
    """
    Test unfolding on the real data
    """
    si_project_dir = si_project_dir(tag)
    folder_name = f'Si_super_deformed{tag}'

    atoms_primitive = read(si_project_dir / f'{folder_name}/Si_prim.cell')
    atoms_supercell = read(si_project_dir / f'{folder_name}/Si_211.cell')
    kpoints, _, labels, _ = read_kpoints(si_project_dir / f'{folder_name}/band.cell', code='castep')

    unfolder: unfold.UnfoldKSet = unfold.UnfoldKSet.from_atoms(np.diag([2, 1, 1]),
                                                               kpoints,
                                                               atoms_primitive,
                                                               atoms_supercell,
                                                               dft_code='castep')
    unfolder.kpoint_labels = labels
    unfolder.write_sc_kpoints(si_project_dir / 'easyunfold_sc_kpoints.cell')

    # Test kpoints generation
    kpoints_sc = read_kpoints(si_project_dir / 'easyunfold_sc_kpoints.cell', code='castep')[0]
    kpoints_sc_ref = read_kpoints(si_project_dir / f'{folder_name}/easyunfold_sc_kpoints.cell', code='castep')[0]
    # shutil.copyfile(si_project_dir / 'easyunfold_sc_kpoints.cell',  # to update
    #                 datapath('Si-project') / f'{folder_name}/easyunfold_sc_kpoints.cell')
    for kpt in kpoints_sc_ref:
        found = False
        for ref in kpoints_sc:
            ref = np.array(ref)
            # Note - allow time reversal symmetry when checking the equivalence of kpoints
            if np.allclose(kpt, ref, atol=1e-9) or np.allclose(kpt, -ref, atol=1e-9):
                found = True
                break
        assert found, f'Kpoint {kpt} not found in the unfolded kpoints'

    # Test unfold
    sws = unfolder.get_spectral_weights(si_project_dir / f'{folder_name}/Si_211_unfold/easyunfold_sc_kpoints.orbitals', ncl=ncl)
    assert sws[0].shape[0] == nspin
    assert len(sws) == len(kpoints)
    assert sws[0].shape[2] == nbands_expected
    assert sws[0].shape[3] == 2

    assert unfolder.is_calculated

    # Spectral weights
    e0, spectral_function = unfolder.get_spectral_function(npoints=500, ncl=ncl)
    assert len(e0) == 500
    assert spectral_function.shape == (nspin, 500, len(kpoints))

    assert unfolder.to_json()


@pytest.mark.parametrize('tag,nspin,ncl, nbands_expected', [('', 1, False, 14)])
def test_unfold_projection(si_project_dir, tag, nspin, ncl, nbands_expected):
    """
    Test unfolding on the real data
    """
    si_project_dir = si_project_dir(tag)
    folder_name = f'Si_super_deformed{tag}'

    atoms_primitive = read(si_project_dir / 'Si/POSCAR')
    atoms_supercell = read(si_project_dir / f'{folder_name}/POSCAR')
    kpoints, _, labels, _ = read_kpoints(si_project_dir / 'KPOINTS_band_low')

    unfolder: unfold.UnfoldKSet = unfold.UnfoldKSet.from_atoms(np.diag([2, 1, 1]), kpoints, atoms_primitive, atoms_supercell)
    unfolder.kpoint_labels = labels
    unfolder.write_sc_kpoints(si_project_dir / 'KPOINTS_sc')

    # Test unfold
    sws = unfolder.get_spectral_weights(si_project_dir / f'{folder_name}/WAVECAR', ncl=ncl)
    assert sws[0].shape[0] == nspin
    assert len(sws) == len(kpoints)
    assert sws[0].shape[2] == nbands_expected
    assert sws[0].shape[3] == 2

    assert unfolder.is_calculated

    unfolder.load_procars(si_project_dir / f'{folder_name}/PROCAR')
    assert 'procars' in unfolder.transient_quantities
    assert 'procars_kmap' in unfolder.transient_quantities


@pytest.mark.parametrize('tag,nspin,ncl,nbands_expected', [
    ('', 1, False, 14),
])
def test_unfold_no_expand(si_project_dir, tag, nspin, ncl, nbands_expected, datapath):
    """
    Test unfolding on the real data without symmetry expansion in the first place
    """
    si_project_dir = si_project_dir(tag)
    folder_name = f'Si_super_deformed{tag}'

    atoms_primitive = read(si_project_dir / 'Si/POSCAR')
    atoms_supercell = read(si_project_dir / f'{folder_name}/POSCAR')
    kpoints, _, labels, _ = read_kpoints(si_project_dir / 'KPOINTS_band_low')

    unfolder: unfold.UnfoldKSet = unfold.UnfoldKSet.from_atoms(np.diag([2, 1, 1]), kpoints, atoms_primitive, atoms_supercell, expand=False)
    unfolder.kpoint_labels = labels
    assert unfolder.nkpts_expand == unfolder.nkpts_orig
    unfolder.write_sc_kpoints(si_project_dir / 'KPOINTS_sc')

    # Test kpoints generation
    kpoints_sc = read_kpoints(si_project_dir / 'KPOINTS_sc')[0]
    # shutil.copyfile(si_project_dir / 'KPOINTS_sc',  # to update
    #                 datapath('Si-project') / f'{folder_name}/KPOINTS_easyunfold')
    kpoints_sc_ref = read_kpoints(si_project_dir / f'{folder_name}/KPOINTS_easyunfold')[0]
    # The kpoints should be a subset of the SC kpoints
    for kpt in kpoints_sc:
        found = False
        for ref in kpoints_sc_ref:
            if np.allclose(kpt, ref):
                found = True
                break
        assert found

    # Test unfold
    sws = unfolder.get_spectral_weights(si_project_dir / f'{folder_name}/WAVECAR', ncl=ncl)
    assert sws[0].shape[0] == nspin
    assert all(x.shape[1] == 1 for x in sws)
    assert len(sws) == len(kpoints)
    assert sws[0].shape[2] == nbands_expected
    assert sws[0].shape[3] == 2

    assert unfolder.is_calculated


def test_atoms_idx_parsing():
    """Test parsing atomic indices"""

    assert unfold.parse_atoms_idx('1,2,3,4') == [0, 1, 2, 3]
    assert unfold.parse_atoms_idx('1-4') == [0, 1, 2, 3]
    assert unfold.parse_atoms_idx('1') == [0]


def test_colourmap():
    """Test handling colour map creation"""
    import matplotlib.pyplot as plt

    cmap = unfold.create_white_colormap('#123456')
    assert cmap(cmap.N) == hex2color('#123456') + (1.0,)

    cmap = unfold.create_white_colormap_from_existing('#123456')
    assert cmap(cmap.N) == hex2color('#123456') + (1.0,)

    cmap = unfold.create_white_colormap_from_existing('Reds')
    assert cmap(cmap.N) == plt.get_cmap('Reds')(256)


def test_atoms_index_process():
    """Test parsing syntax for passing atomic indices"""

    assert unfold.parse_atoms_idx('1,2') == [0, 1]
    assert unfold.parse_atoms_idx('1, 2') == [0, 1]
    assert unfold.parse_atoms_idx('1 , 2') == [0, 1]
    assert unfold.parse_atoms_idx('1,2,4-5') == [0, 1, 3, 4]
    assert unfold.parse_atoms_idx('1-10') == list(range(0, 10))


def test_project_option_processing():
    """Test for parsing projection options """

    a, b = unfold.process_projection_options('1,2', None)
    assert a == [0, 1]
    assert b == 'all'

    a, b = unfold.process_projection_options('1,2', 's,p ')
    assert a == [0, 1]
    assert b == ['s', 'p']
