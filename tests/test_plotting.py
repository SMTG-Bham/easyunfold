"""
Test the effective mass functionality
"""
from pathlib import Path
import pytest
import numpy as np
from monty.serialization import loadfn
from ase.io import read
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from easyunfold.unfold import UnfoldKSet
from easyunfold.utils import read_kpoints
from easyunfold.effective_mass import EffectiveMass
import easyunfold.plotting as pl

# pylint:disable=redefined-outer-name


@pytest.fixture
def silicon_unfolded(si_project_dir) -> UnfoldKSet:
    """
    Test unfolding on the real data
    """
    tag = ''
    ncl = False
    si_project_dir = si_project_dir(tag)
    folder_name = f'Si_super_deformed{tag}'

    atoms_primitive = read(si_project_dir / 'Si/POSCAR')
    atoms_supercell = read(si_project_dir / f'{folder_name}/POSCAR')
    kpoints, _, labels, _ = read_kpoints(si_project_dir / 'KPOINTS_band_low')

    unfolder: UnfoldKSet = UnfoldKSet.from_atoms(np.diag([2, 1, 1]), kpoints, atoms_primitive, atoms_supercell)
    unfolder.kpoint_labels = labels
    # Test unfold
    unfolder.get_spectral_weights(si_project_dir / f'{folder_name}/WAVECAR', ncl=ncl)
    return unfolder


@pytest.fixture(scope='module')
def unfold_obj() -> UnfoldKSet:
    """Return an unfolding object"""
    return loadfn(Path(__file__).parent / 'test_data/mgo.json')


def test_plotting(unfold_obj: UnfoldKSet):
    """Test the plotting routine"""

    plotter = pl.UnfoldPlotter(unfold_obj)
    fig = plotter.plot_spectral_weights()
    assert isinstance(fig, Figure)

    eng, sf = unfold_obj.get_spectral_function(npoints=200)
    fig = plotter.plot_spectral_function(eng, sf)
    assert isinstance(fig, Figure)

    # Plot onto existing figures
    fig, ax = plt.subplots(1, 1)
    fig = plotter.plot_spectral_weights(ax=ax)
    assert isinstance(fig, Figure)

    eng, sf = unfold_obj.get_spectral_function(npoints=200)
    fig = plotter.plot_spectral_function(eng, sf)
    assert isinstance(fig, Figure)


def test_plotting_projection(unfold_obj: UnfoldKSet):
    """Test plotting with projections from PROCAR"""
    # Test using PROCAR for projection plot
    plotter = pl.UnfoldPlotter(unfold_obj)
    procar_path = Path(__file__).parent / 'test_data/PROCAR.mgo'
    fig = plotter.plot_projected(procar_path, atoms_idx='0,1|2,3', npoints=200, use_subplot=False)
    assert isinstance(fig, Figure)

    fig = plotter.plot_projected(procar_path, atoms_idx='0,1|2,3', npoints=200, use_subplot=True)
    assert isinstance(fig, Figure)

    # test options
    poscar_path = Path(__file__).parent / 'test_data/POSCAR.mgo'
    fig = plotter.plot_projected(procar_path,
                                 eref=2,
                                 gamma=False,
                                 ncl=False,
                                 npoints=200,
                                 sigma=0.15,
                                 symm_average=False,
                                 figsize=(2, 3),
                                 ylim=(-2, 5),
                                 dpi=500,
                                 vscale=2.0,
                                 contour_plot=True,
                                 alpha=0.4,
                                 save='test.png',
                                 ax=None,
                                 cmap='PuBu',
                                 show=True,
                                 title='Test',
                                 atoms='Mg,O',
                                 poscar=poscar_path,
                                 orbitals='s,p',
                                 use_subplot=True,
                                 colours='r,g',
                                 colorspace='luvlch',
                                 intensity=0.5)
    assert isinstance(fig, Figure)


def test_color_interpolation():
    """Test interpolating colours"""
    colours = ['r', 'g', 'b', 'cyan']
    weights = np.random.rand(10, 4)
    output = pl.interpolate_colors(colours, weights)

    assert output.shape == (10, 3)
    assert np.max(output) <= 1.0


def test_effective_mass_plot(silicon_unfolded):
    """Test plotting effective mass information"""

    unfold: UnfoldKSet = silicon_unfolded
    plotter = pl.UnfoldPlotter(unfold)
    eff = EffectiveMass(unfold)
    engs, sf = unfold.get_spectral_function(npoints=200)
    fig = plotter.plot_effective_mass(eff, engs, sf)

    assert isinstance(fig, Figure)
