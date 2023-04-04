"""
Compatibility layers for handling wavefunctions
"""

# pylint: disable= protected-access, useless-super-delegation
import numpy as np
from ase.units import Hartree

from castepxbin.wave import WaveFunction as CastepWF
from .wavecar import Wavecar


class WaveFunction:
    """
    Interface for accessing wavefunction data

    All indexings are one-base rather than zero-based as in python
    """

    def __init__(self, wfc):
        self.wfc = wfc

    @property
    def kpoints(self):
        """Kpoints as row vectors"""
        raise NotImplementedError

    @property
    def nkpts(self):
        return self.kpoints.shape[0]

    @property
    def nspins(self):
        """Number of spins"""
        raise NotImplementedError

    @property
    def mesh_size(self):
        raise NotImplementedError

    @property
    def bands(self):
        """KS band energies in shape (ns, nk, nb)"""
        raise NotImplementedError

    @property
    def nbands(self):
        """Number of KS bands"""
        return self.bands.shape[-1]

    @property
    def occupancies(self):
        """Occupancies of each band"""
        raise NotImplementedError

    def get_gvectors(self, ik):
        """Return the gvectors at a kpoint with shape (nwaves, 3)"""
        raise NotImplementedError

    def get_band_coeffs(self, ispin, ik, ib, norm=True):
        """Return band coefficients"""
        raise NotImplementedError


class VaspWaveFunction(WaveFunction):
    """Interface for accessing WAVECAR"""

    def __init__(self, wfc: Wavecar):
        super().__init__(wfc)

    @property
    def kpoints(self):
        return self.wfc._kvecs

    @property
    def occupancies(self):
        return self.wfc._occs

    @property
    def nspins(self):
        return self.wfc._nspin

    @property
    def mesh_size(self):
        return self.wfc._ngrid

    @property
    def bands(self):
        """KS band energies in shape (ns, nk, nb)"""
        return self.wfc._bands

    def get_gvectors(self, ik):
        """Return the gvectors at a kpoint"""
        return self.wfc.get_gvectors(ik)

    def get_band_coeffs(self, ispin, ik, ib, norm=True):
        """Return plane wave coefficients at a specific band"""
        return self.wfc.read_band_coeffs(ispin, ik, ib, norm=norm)


class CastepWaveFunction(WaveFunction):
    """
    Interface for reading wave function data from a CASTEP calculation.
    """

    def __init__(self, wfc: CastepWF):
        super().__init__(wfc)
        self.wfc: CastepWF

    @classmethod
    def from_file(cls, fname):
        return cls(CastepWF.from_file(fname))

    @property
    def kpoints(self):
        return self.wfc.kpts.T

    @property
    def occupancies(self):
        return self.wfc.occupancies.T

    @property
    def mesh_size(self):
        return self.wfc.mesh_size

    @property
    def bands(self):
        return self.wfc.eigenvalues.T * Hartree

    @property
    def nspins(self):
        return self.wfc.nspins

    def get_gvectors(self, ik):
        """Return the G-vector at a kpoint"""
        return self.wfc.get_gvectors(ik - 1).T

    def get_band_coeffs(self, ispin, ik, ib, norm=True):
        """Return the plane wave coefficients for a band"""
        coeffs = self.wfc.get_plane_wave_coeffs(ispin - 1, ik - 1, ib - 1)
        if norm:
            coeffs = coeffs / np.linalg.norm(coeffs)
        return coeffs
