"""
Code for reading the eigenvalues and plane wave coefficients from the WAVECAR file.
"""

import os
from math import sqrt

import numpy as np
from .vasp_constant import RYTOEV, TPI, AUTOA, HSQDTM


class Wavecar:  # pylint: disable=too-many-instance-attributes
    """
    Class for processing VASP Pseudowavefunction stored in WAVECAR.

    This class is a trimmed-down from that of [PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding)
    by QijingZheng <zqj.kaka@gmail.com>.

    Only the functions needed for band unfolding remains.
    The original class has many other functionalities such as ELF calculation.

    The format of VASP WAVECAR, as shown in
        http://www.andrew.cmu.edu/user/feenstra/wavetrans/
    is:
    ```
        Record-length #spin components RTAG(a value specifying the precision)
        #k-points #bands ENCUT(maximum energy for plane waves)
        LatVec-A
        LatVec-B
        LatVec-C
        Loop over spin
           Loop over k-points
              #plane waves, k vector
              Loop over bands
                 band energy, band occupation
              End loop over bands
              Loop over bands
                 Loop over plane waves
                    Plane-wave coefficient
                 End loop over plane waves
              End loop over bands
           End loop over k-points
        End loop over spin
    ```
    """

    def __init__(self, fnm='WAVECAR', lsorbit=False, lgamma=False, gamma_half='x'):
        r"""
        Initialization.

        :param fnm: Path to the `WAVECAR` file
        :param lsorbit: Whether the `WAVECAR` is from a spint-orbit coupling calculation
        :param lgamma: Wether the `WAVECAR` is from a $\Gamma$-only calculation
        """

        self._fname = fnm
        # the directory containing the input file
        self._dname = os.path.dirname(fnm)
        if self._dname == '':
            self._dname = '.'

        self._lsoc = lsorbit
        self._lgam = lgamma
        self._gam_half = gamma_half.lower()

        assert not (lsorbit and lgamma), 'Cannot have both `lsorbit` and `lgamma`!'
        assert self._gam_half in ['x', 'z'], 'Gamma_half must be "x" or "z"'

        self._wfc = open(self._fname, 'rb')  # pylint:disable=consider-using-with

        # read the basic information
        self.read_wf_header()
        # read the band information
        self.read_bands()

        if self._lsoc:
            assert self._nspin == 1, 'NSPIN = 1 for noncollinear version WAVECAR!'

    def close(self):
        """Close the file handle"""
        self._wfc.close()

    def is_soc(self):
        """
        Return `True` if the WAVECAR from an SOC calculation.
        """
        return bool(self._lsoc)

    def is_gamma(self):
        r"""
        Return `True` is the WAVECAR is from an $\Gamma$-only calculation.
        """
        return bool(self._lgam)

    def read_wf_header(self):
        """
        Read the system information from WAVECAR, which is written in the first
        two record.

        :rec1:
            `recl, nspin, rtag`
        :rec2:
            `nkpts, nbands, encut, ((cell(i,j) i=1, 3), j=1, 3)`
        """

        # goto the start of the file and read the first record
        self._wfc.seek(0)
        self._recl, self._nspin, self._rtag = np.array(np.fromfile(self._wfc, dtype=np.float64, count=3), dtype=np.int64)
        self._precision = self.set_wf_prec()
        # the second record
        self._wfc.seek(self._recl)
        dump = np.fromfile(self._wfc, dtype=np.float64, count=12)

        self._nkpts = int(dump[0])  # No. of k-points
        self._nbands = int(dump[1])  # No. of bands
        self._encut = dump[2]  # Energy cutoff
        # real space supercell basis
        self._realspace_cell = dump[3:].reshape((3, 3))
        # real space supercell volume
        self._realspace_cell_volume = np.linalg.det(self._realspace_cell)
        # reciprocal space supercell volume
        self._reciprocal_cell_volume = np.linalg.inv(self._realspace_cell).T

        # Minimum FFT grid size
        cell_abc = np.linalg.norm(self._realspace_cell, axis=1)
        cut_off = np.ceil(sqrt(self._encut / RYTOEV) / (TPI / (cell_abc / AUTOA)))
        self._ngrid = np.array(2 * cut_off + 1, dtype=int)

    def set_wf_prec(self):
        """
        Set wavefunction coefficients precision:

        :TAG = 45200:
            single precision complex, np.complex64, or complex(qs)
        :TAG = 45210:
            double precision complex, np.complex128, or complex(q)
        """
        if self._rtag == 45200:
            return np.complex64
        if self._rtag == 45210:
            return np.complex128
        if self._rtag == 53300:
            raise ValueError('VASP5 WAVECAR format, not implemented yet')
        if self._rtag == 53310:
            raise ValueError('VASP5 WAVECAR format with double precision ' + 'coefficients, not implemented yet')
        raise ValueError(f'Invalid TAG values: {self._rtag}')

    def read_bands(self) -> tuple:
        """
        Extract KS energies and Fermi occupations from WAVECAR.
        """

        self._nplws = np.zeros(self._nkpts, dtype=int)
        self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
        self._bands = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)
        self._occs = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)

        for ii in range(self._nspin):
            for jj in range(self._nkpts):
                rec = self.locate_rec(ii + 1, jj + 1, 1) - 1
                self._wfc.seek(rec * self._recl)
                dump = np.fromfile(self._wfc, dtype=np.float64, count=4 + 3 * self._nbands)
                if ii == 0:
                    self._nplws[jj] = int(dump[0])
                    self._kvecs[jj] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self._bands[ii, jj, :] = dump[:, 0]
                self._occs[ii, jj, :] = dump[:, 2]

        if self._nkpts > 1:
            tmp = np.linalg.norm(np.dot(np.diff(self._kvecs, axis=0), self._reciprocal_cell_volume), axis=1)
            self._kpath = np.concatenate(([
                0,
            ], np.cumsum(tmp)))
        else:
            self._kpath = None
        return self._kpath, self._bands

    def get_gvectors(self, ikpt=1, force_gamma=False, check_consistency=True):
        r"""
        Generate the G-vectors that satisfies the following relation

        $$
            \frac{|\vec{G} + \vec{k}|^2}{2} < E_{cut}
        $$

        """
        assert 1 <= ikpt <= self._nkpts, 'Invalid kpoint index!'

        kvec = self._kvecs[ikpt - 1]
        # force_gamma: consider gamma-only case regardless of the actual setting
        lgam = True if force_gamma else self._lgam

        fx, fy, fz = [np.arange(n, dtype=int) for n in self._ngrid]
        fx[self._ngrid[0] // 2 + 1:] -= self._ngrid[0]
        fy[self._ngrid[1] // 2 + 1:] -= self._ngrid[1]
        fz[self._ngrid[2] // 2 + 1:] -= self._ngrid[2]
        if lgam:
            if self._gam_half == 'x':
                fx = fx[:self._ngrid[0] // 2 + 1]
            else:
                fz = fz[:self._ngrid[2] // 2 + 1]

    # In meshgrid, fx run the fastest, fz the slowest
        gz, gy, gx = np.array(np.meshgrid(fz, fy, fx, indexing='ij')).reshape((3, -1))
        kgrid = np.array([gx, gy, gz], dtype=float).T
        if lgam:
            if self._gam_half == 'z':
                kgrid = kgrid[(gz > 0) | ((gz == 0) & (gy > 0)) | ((gz == 0) & (gy == 0) & (gx >= 0))]
            else:
                kgrid = kgrid[(gx > 0) | ((gx == 0) & (gy > 0)) | ((gx == 0) & (gy == 0) & (gz >= 0))]

        # Kinetic_Energy = (G + k)**2 / 2
        kinetic_energy = HSQDTM * np.linalg.norm(np.dot(kgrid + kvec[np.newaxis, :], TPI * self._reciprocal_cell_volume), axis=1)**2
        # find Gvectors where (G + k)**2 / 2 < ENCUT
        g_vectors = kgrid[np.where(kinetic_energy < self._encut)[0]]

        # Check if the calculated number of planewaves and the one recorded in the
        # WAVECAR are equal
        if check_consistency:

            if g_vectors.shape[0] != self._nplws[ikpt - 1]:
                if g_vectors.shape[0] * 2 == self._nplws[ikpt - 1]:
                    if not self._lsoc:
                        raise ValueError("""
                        It seems that you are reading a WAVECAR from a NONCOLLINEAR VASP.
                        Please set 'lsorbit = True' when loading the WAVECAR.
                        For example:

                            wfc = Wavecar('WAVECAR', lsorbit=True)
                        """)
                elif g_vectors.shape[0] == 2 * self._nplws[ikpt - 1] - 1:
                    if not self._lgam:
                        raise ValueError("""
                        It seems that you are reading a WAVECAR from a GAMMA-ONLY VASP.  Please set
                        'lgamma = True' when loading the WAVECAR.  Moreover, you may want to set
                        "gamma_half" if you are using VASP version <= 5.2.x.  For VASP <= 5.2.x, check
                        which FFT VASP uses by the following command:

                            $ grep 'use.* FFT for wave' OUTCAR

                        Then

                            # for parallel FFT, VASP <= 5.2.x
                            wfc = Wavecar('WAVECAR', lgamma=True, gamma_half='z')

                            # for serial FFT, VASP <= 5.2.x
                            wfc = Wavecar('WAVECAR', lgamma=True, gamma_half='x')

                        For VASP >= 5.4, WAVECAR is written with x-direction half grid regardless of
                        parallel or serial FFT.

                            # "gamma_half" default to "x" for VASP >= 5.4
                            wfc = Wavecar('WAVECAR', lgamma=True, gamma_half='x')
                        """)
                else:
                    raise ValueError(f"""
                    NO. OF PLANEWAVES NOT CONSISTENT:

                        THIS CODE -> {g_vectors.shape[0]}
                        FROM VASP -> {self._nplws[ikpt - 1] // 2 if self._lsoc else self._nplws[ikpt - 1]}
                           NGRIDS -> {np.prod(self._ngrid)}
                    """)

        return np.asarray(g_vectors, dtype=int)

    def read_band_coeffs(self, ispin=1, ikpt=1, iband=1, norm=False):
        """
        Read the planewave coefficients of specified KS states.
        """

        self.check_index(ispin, ikpt, iband)

        rec = self.locate_rec(ispin, ikpt, iband)
        self._wfc.seek(rec * self._recl)

        nplw = self._nplws[ikpt - 1]
        dump = np.fromfile(self._wfc, dtype=self._precision, count=nplw)

        cg = np.asarray(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def locate_rec(self, ispin=1, ikpt=1, iband=1):
        """
        Return the rec position for specified KS state.
        """

        self.check_index(ispin, ikpt, iband)

        rec = 2 + (ispin - 1) * self._nkpts * (self._nbands + 1) + \
                  (ikpt - 1) * (self._nbands + 1) + \
            iband
        return rec

    def check_index(self, ispin, ikpt, iband):
        """
        Check if the index is valid!
        """
        assert 1 <= ispin <= self._nspin, 'Invalid spin index!'
        assert 1 <= ikpt <= self._nkpts, 'Invalid kpoint index!'
        assert 1 <= iband <= self._nbands, 'Invalid band index!'
