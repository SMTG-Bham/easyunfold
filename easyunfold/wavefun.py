"""
Compatibility layers for handling wavefunctions
"""

# pylint: disable= protected-access


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
