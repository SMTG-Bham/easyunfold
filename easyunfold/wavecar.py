"""
Code for reading the eigenvalues and plane wave coefficients from the WAVECAR file.
"""

import os
from math import sqrt

import numpy as np
from .vasp_constant import RYTOEV, TPI, AUTOA, HSQDTM


class Wavecar:
    """
    Class for processing VASP Pseudowavefunction stored in WAVECAR.

    This class is a trimmed-down from that of [PyVaspwfc](https://github.com/QijingZheng/VaspBandUnfolding)
    by QijingZheng <zqj.kaka@gmail.com>.

    Only the functions needed for band unfolding remains.
    The original class has many other functionalities such as ELF calculation.

    The format of VASP WAVECAR, as shown in
        http://www.andrew.cmu.edu/user/feenstra/wavetrans/
    is:
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
    """

    def __init__(self, fnm='WAVECAR', lsorbit=False, lgamma=False, gamma_half='x', omp_num_threads=1):
        """
        Initialization.
        """

        self._fname = fnm
        # the directory containing the input file
        self._dname = os.path.dirname(fnm)
        if self._dname == '':
            self._dname = '.'

        self._lsoc = lsorbit
        self._lgam = lgamma
        self._gam_half = gamma_half.lower()

        # It seems that some modules in scipy uses OPENMP, it is therefore
        # desirable to set the OMP_NUM_THREADS to tune the parallization.
        os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)

        assert not (lsorbit and lgamma), 'The two settings conflict!'
        assert self._gam_half == 'x' or self._gam_half == 'z', \
            'Gamma_half must be "x" or "z"'

        try:
            self._wfc = open(self._fname, 'rb')
        except:
            raise IOError('Failed to open %s' % self._fname)

        # read the basic information
        self.readWFHeader()
        # read the band information
        self.readWFBand()

        if self._lsoc:
            assert self._nspin == 1, 'NSPIN = 1 for noncollinear version WAVECAR!'

    def isSocWfc(self):
        """
        Is the WAVECAR from an SOC calculation?
        """
        return True if self._lsoc else False

    def isGammaWfc(self):
        """
        Is the WAVECAR from an SOC calculation?
        """
        return True if self._lgam else False

    def readWFHeader(self):
        """
        Read the system information from WAVECAR, which is written in the first
        two record.

        rec1: recl, nspin, rtag
        rec2: nkpts, nbands, encut, ((cell(i,j) i=1, 3), j=1, 3)
        """

        # goto the start of the file and read the first record
        self._wfc.seek(0)
        self._recl, self._nspin, self._rtag = np.array(np.fromfile(self._wfc, dtype=np.float64, count=3), dtype=np.int64)
        self._WFPrec = self.setWFPrec()
        # the second record
        self._wfc.seek(self._recl)
        dump = np.fromfile(self._wfc, dtype=np.float64, count=12)

        self._nkpts = int(dump[0])  # No. of k-points
        self._nbands = int(dump[1])  # No. of bands
        self._encut = dump[2]  # Energy cutoff
        # real space supercell basis
        self._Acell = dump[3:].reshape((3, 3))
        # real space supercell volume
        self._Omega = np.linalg.det(self._Acell)
        # reciprocal space supercell volume
        self._Bcell = np.linalg.inv(self._Acell).T

        # Minimum FFT grid size
        Anorm = np.linalg.norm(self._Acell, axis=1)
        CUTOF = np.ceil(sqrt(self._encut / RYTOEV) / (TPI / (Anorm / AUTOA)))
        self._ngrid = np.array(2 * CUTOF + 1, dtype=int)

    def setWFPrec(self):
        """
        Set wavefunction coefficients precision:
            TAG = 45200: single precision complex, np.complex64, or complex(qs)
            TAG = 45210: double precision complex, np.complex128, or complex(q)
        """
        if self._rtag == 45200:
            return np.complex64
        elif self._rtag == 45210:
            return np.complex128
        elif self._rtag == 53300:
            raise ValueError('VASP5 WAVECAR format, not implemented yet')
        elif self._rtag == 53310:
            raise ValueError('VASP5 WAVECAR format with double precision ' + 'coefficients, not implemented yet')
        else:
            raise ValueError('Invalid TAG values: {}'.format(self._rtag))

    def readWFBand(self):
        """
        Extract KS energies and Fermi occupations from WAVECAR.
        """

        self._nplws = np.zeros(self._nkpts, dtype=int)
        self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
        self._bands = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)
        self._occs = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)

        for ii in range(self._nspin):
            for jj in range(self._nkpts):
                rec = self.whereRec(ii + 1, jj + 1, 1) - 1
                self._wfc.seek(rec * self._recl)
                dump = np.fromfile(self._wfc, dtype=np.float64, count=4 + 3 * self._nbands)
                if ii == 0:
                    self._nplws[jj] = int(dump[0])
                    self._kvecs[jj] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self._bands[ii, jj, :] = dump[:, 0]
                self._occs[ii, jj, :] = dump[:, 2]

        if self._nkpts > 1:
            tmp = np.linalg.norm(np.dot(np.diff(self._kvecs, axis=0), self._Bcell), axis=1)
            self._kpath = np.concatenate(([
                0,
            ], np.cumsum(tmp)))
        else:
            self._kpath = None
        return self._kpath, self._bands

    def gvectors(self, ikpt=1, force_Gamma=False, check_consistency=True):
        """
        Generate the G-vectors that satisfies the following relation
            (G + k)**2 / 2 < ENCUT
        """
        assert 1 <= ikpt <= self._nkpts, 'Invalid kpoint index!'

        kvec = self._kvecs[ikpt - 1]
        # force_Gamma: consider gamma-only case regardless of the actual setting
        lgam = True if force_Gamma else self._lgam

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
        KENERGY = HSQDTM * np.linalg.norm(np.dot(kgrid + kvec[np.newaxis, :], TPI * self._Bcell), axis=1)**2
        # find Gvectors where (G + k)**2 / 2 < ENCUT
        Gvec = kgrid[np.where(KENERGY < self._encut)[0]]

        # Check if the calculated number of planewaves and the one recorded in the
        # WAVECAR are equal
        if check_consistency:

            if Gvec.shape[0] != self._nplws[ikpt - 1]:
                if Gvec.shape[0] * 2 == self._nplws[ikpt - 1]:
                    if not self._lsoc:
                        raise ValueError("""
                        It seems that you are reading a WAVECAR from a NONCOLLINEAR VASP.
                        Please set 'lsorbit = True' when loading the WAVECAR.
                        For example:

                            wfc = Wavecar('WAVECAR', lsorbit=True)
                        """)
                elif Gvec.shape[0] == 2 * self._nplws[ikpt - 1] - 1:
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
                    raise ValueError("""
                    NO. OF PLANEWAVES NOT CONSISTENT:

                        THIS CODE -> %d
                        FROM VASP -> %d
                           NGRIDS -> %d
                    """ % (Gvec.shape[0], self._nplws[ikpt - 1] // 2 if self._lsoc else self._nplws[ikpt - 1], np.prod(self._ngrid)))

        return np.asarray(Gvec, dtype=int)

    def readBandCoeff(self, ispin=1, ikpt=1, iband=1, norm=False):
        """
        Read the planewave coefficients of specified KS states.
        """

        self.checkIndex(ispin, ikpt, iband)

        rec = self.whereRec(ispin, ikpt, iband)
        self._wfc.seek(rec * self._recl)

        nplw = self._nplws[ikpt - 1]
        dump = np.fromfile(self._wfc, dtype=self._WFPrec, count=nplw)

        cg = np.asarray(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def whereRec(self, ispin=1, ikpt=1, iband=1):
        """
        Return the rec position for specified KS state.
        """

        self.checkIndex(ispin, ikpt, iband)

        rec = 2 + (ispin - 1) * self._nkpts * (self._nbands + 1) + \
                  (ikpt - 1) * (self._nbands + 1) + \
            iband
        return rec

    def checkIndex(self, ispin, ikpt, iband):
        """
        Check if the index is valid!
        """
        assert 1 <= ispin <= self._nspin, 'Invalid spin index!'
        assert 1 <= ikpt <= self._nkpts, 'Invalid kpoint index!'
        assert 1 <= iband <= self._nbands, 'Invalid band index!'
