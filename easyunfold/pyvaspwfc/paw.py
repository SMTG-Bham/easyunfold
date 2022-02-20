#!/usr/bin/env python

import re
import numpy as np
from .vasp_constant import *
from .sph_harm import sph_r, sph_c


def fftchk1(n):
    """
    Check if n can be factorized into products of 2, 3 and 5.

    From VASP fft3dcray.F, FFTCHK1
    """
    if n % 2 == 1:
        return False

    nmax = np.array(
        np.log(n) / np.log([2, 3, 5]),
        dtype=int
    )
    ijk = np.mgrid[
        0:nmax[2]+1,
        0:nmax[1]+1,
        1:nmax[0]+1,     # n should be even number
    ].reshape((3, -1)).T

    for i, j, k in ijk:
        if (5**i * 3**j * 2**k) == n:
            return True
    return False


def fftchk(ngrid):
    '''
    Return the next correct settings for 3d FFT.

    From VASP fft3dcray.F, FFTCHK
    '''

    ngrid = np.asarray(ngrid, dtype=int)
    assert ngrid.shape == (3,)

    for ii in range(3):
        while not fftchk1(ngrid[ii]):
            ngrid[ii] += 1

    return ngrid


def gvectors(cell, encut, kvec, ngrid=None,
             lgam=False, gamma_half='x', force_gamma=False,
             lsoc=False):
    '''
    Generate the G-vectors that satisfies the following relation
        (G + k)**2 / 2 < ENCUT
    '''
    # Minimum FFT grid size
    Bcell = np.linalg.inv(cell).T     # reciprocal space supercell volume
    if ngrid is None:
        Anorm = np.linalg.norm(cell, axis=1)
        CUTOF = np.ceil(
            np.sqrt(encut / RYTOEV) / (TPI / (Anorm / AUTOA))
        )
        ngrid = np.array(2 * CUTOF + 1, dtype=int)

    kvec = np.asarray(kvec)

    # force_Gamma: consider gamma-only case regardless of the actual setting
    if force_gamma:
        lgam = True

    fx, fy, fz = [np.arange(n, dtype=int) for n in ngrid]
    fx[ngrid[0] // 2 + 1:] -= ngrid[0]
    fy[ngrid[1] // 2 + 1:] -= ngrid[1]
    fz[ngrid[2] // 2 + 1:] -= ngrid[2]
    if lgam:
        if gam_half == 'x':
            fx = fx[:ngrid[0] // 2 + 1]
        else:
            fz = fz[:ngrid[2] // 2 + 1]

    gz, gy, gx = np.array(
        np.meshgrid(fz, fy, fx, indexing='ij')
    ).reshape((3, -1))
    kgrid = np.array([gx, gy, gz], dtype=float).T
    if lgam:
        if gamma_half == 'z':
            kgrid = kgrid[
                (gz > 0) |
                ((gz == 0) & (gy > 0)) |
                ((gz == 0) & (gy == 0) & (gx >= 0))
            ]
        else:
            kgrid = kgrid[
                (gx > 0) |
                ((gx == 0) & (gy > 0)) |
                ((gx == 0) & (gy == 0) & (gz >= 0))
            ]

    # Kinetic_Energy = (G + k)**2 / 2
    # HSQDTM    =  hbar**2/(2*ELECTRON MASS)
    KENERGY = HSQDTM * np.linalg.norm(
        np.dot(kgrid + kvec[np.newaxis, :], TPI*Bcell), axis=1
    )**2
    # find Gvectors where (G + k)**2 / 2 < ENCUT
    Gvec = kgrid[np.where(KENERGY < encut)[0]]

    return np.asarray(Gvec, dtype=int)


class pawpotcar(object):
    '''
    Read projector functions and ae/ps partialwaves from VASP PBE POTCAR.
    '''

    NPSQNL = 100      # no. of data for projectors in reciprocal space
    NPSRNL = 100      # no. of data for projectors in real space

    def __init__(self, potstr=None, potfile=None):
        '''
        PAW POTCAR provides the PAW projector functions in real and reciprocal
        space. 
        '''
        if (potfile is not None) and (potstr is None):
            # read POTCAR of the first element if there are more than one in the
            # potfile
            potstr = open(potfile).read().split('End of Dataset')[0]
        assert len(potstr.strip()) != 0, "POTCAR string should not be empty!"

        non_radial_part, radial_part = potstr.split('PAW radial sets', 1)

        # read the projector functions in real/reciprocal space
        self.read_proj(non_radial_part)
        # read the ae/ps partial waves in the core region
        self.read_partial_wfc(radial_part)
        # c-spline interpolation of the projector function
        self.csplines()

    def read_proj(self, datastr):
        '''
        Read the projector functions in reciprocal space.
        '''
        dump = datastr.split('Non local Part')
        head = dump[0].strip().split('\n')
        non_local_part = dump[1:]

        # element of the potcar
        self.element = head[0].split()[1]
        # maximal G for reciprocal non local projectors
        self.proj_gmax = float(head[-1].split()[0])

        qprojs = []
        rprojs = []
        proj_l = []
        for proj in non_local_part:
            dump = proj.split('Reciprocal Space Part')
            ln_rmax_dion_part = dump[0].strip().split()
            l, nlproj = [int(xx) for xx in ln_rmax_dion_part[:2]]

            proj_l += [l] * nlproj
            self.proj_rmax = float(ln_rmax_dion_part[2])

            dion = np.asarray(ln_rmax_dion_part[3:], dtype=float)

            for rr in dump[1:]:
                reci, real = rr.split('Real Space Part')
                qprojs.append(
                    np.fromstring(reci, float, sep=' ')
                )
                rprojs.append(
                    np.fromstring(real, float, sep=' ')
                )

        # the real space radial grid for the projector functions
        self.proj_rgrid = np.arange(self.NPSRNL) * self.proj_rmax / self.NPSRNL
        # the reciprocal space radial grid for the projector functions
        self.proj_qgrid = np.arange(self.NPSQNL) * self.proj_gmax / self.NPSQNL
        # L quantum number for each projector functions
        self.proj_l = np.asarray(proj_l, dtype=int)
        # projector functions in reciprocal space
        self.qprojs = np.asarray(qprojs, dtype=float)
        # projector functions in real space
        self.rprojs = np.asarray(rprojs, dtype=float)

    def read_partial_wfc(self, datastr):
        '''
        Read the ps/ae partial waves.

        The data structure in POTCAR:

             grid
             aepotential
             core charge-density
             kinetic energy-density
             mkinetic energy-density pseudized
             local pseudopotential core
             pspotential valence only
             core charge-density (pseudized)
             pseudo wavefunction
             ae wavefunction
             ...
             pseudo wavefunction
             ae wavefunction
        '''
        data = datastr.strip().split('\n')
        nmax = int(data[0].split()[0])
        grid_start_idx = data.index(" grid") + 1

        core_data = np.array([
            x for line in data[grid_start_idx:]
            for x in line.strip().split()
            if not re.match(r'\ \w+', line)
        ], dtype=float)
        core_data = core_data.reshape((-1, nmax))
        # number of projectors
        nproj = self.proj_l.size

        # core region logarithmic radial grid
        self.rgrid = core_data[0]
        # core region all-electron potential
        self.paw_aepot = core_data[1]
        # core region pseudo wavefunction
        self.paw_ps_wfc = core_data[-nproj*2::2, :]
        # core region all-electron wavefunctions
        self.paw_ae_wfc = core_data[-nproj*2+1::2, :]

    def csplines(self):
        '''
        Cubic spline interpolation of both the real and reciprocal space
        projector functions.
        '''

        from scipy.interpolate import CubicSpline as cs

        # for reciprocal space projector functions, natural boundary condition
        # (Y'' = 0) is applied at both ends.
        self.spl_qproj = [
            cs(self.proj_qgrid, qproj, bc_type='natural') for qproj in
            self.qprojs
        ]
        # for real space projector functions, natural boundary condition
        # (Y'' = 0) is applied at the point N.
        self.spl_rproj = []
        for l, rproj in zip(self.proj_l, self.rprojs):
            # Copy from VASP pseudo.F ln. 444 - 448, I don't know why y1p depend on "l".
            if l == 1:
                y1p = (rproj[1] - rproj[0]) / (self.proj_rmax / self.NPSRNL)
            else:
                y1p = 0.0
            self.spl_rproj.append(
                cs(self.proj_rgrid, rproj, bc_type=((1, y1p), (2, 0)))
            )

    def set_simpi_weight(self):
        '''
        Setup weights for simpson integration on radial grid any radial integral
        can then be evaluated by just summing all radial grid points with the
        weights SI

        \int dr = \sum_i w(i) * f(i)
        '''

        self.rad_simp_w = np.zeros_like(self.rgrid)
        # Number of points in the radial grid
        N = self.rgrid.size
        # Logarithmic grid: R(i+1) / R(i) = exp(H)
        # Logarithmic grid: R(i) = R(0) * exp(H*i)
        H = np.log((self.rgrid[-1] / self.rgrid[0])**(1. / (N-1)))

        for ii in range(N-1, 1, -2):
            self.rad_simp_w[ii] = H * self.rgrid[ii] / \
                3. + self.rad_simp_w[ii]
            self.rad_simp_w[ii-1] = H * self.rgrid[ii-1] * 4. / 3.
            self.rad_simp_w[ii-2] = H * self.rgrid[ii-2] / 3.

    def radial_simp_int(self, f):
        '''
        Simpson integration of a function on the logarithmic radial grid.
        '''
        if not hasattr(self, "rad_simp_w"):
            self.set_simpi_weight()
        f = np.asarray(f)

        return np.sum(self.rad_simp_w * f)

    def get_Qij(self):
        '''
        Calculate the quantity

            Q_{ij} = < \phi_i^{AE} | \phi_j^{AE} > - < phi_i^{PS} | phi_j^{PS} >

        where \phi^{AE/PS} are the PAW AE/PS waves, which are real functions in
        VASP PAW POTCAR.

        In POTCAR, only the radial part of the AE/PS partial waves are stored.
        In order to get the total AE/PS waves, a spherical harmonics should be
        multiplied to the radial part, i.e. 

            \psi^{AE/PS}(r) = (1. / r) * \phi^{AE/PS}(r) * Y(l, m)

        where \psi(r) is the total AE/PS partial waves, and \phi(r) are the ones
        stored in POTCAR. Y(l, m) is the real spherical harmonics with quantum
        number l and m. Note the "1 / r" in front of the equation. In practice,
        the r**(-2) term disappears when multiply with the integration volume
        "r**2 sin(theta) d theta d phi".

        In theory, the "Qij" should be integrated inside the PAW cut-off radius,
        but since the AE and PS partial waves are indentical outside the cut-off
        radius, therefore the terms outside the cut-off radius cancel each
        other. 
        '''

        if not hasattr(self, 'paw_qij'):
            # The l quantum number of each partial waves
            LL = [l for l in self.proj_l for m in range(-l, l+1)]
            # The m quantum number of each partial waves
            MM = [m for l in self.proj_l for m in range(-l, l+1)]
            # The index of each partial waves
            WW = [i for i, l in enumerate(self.proj_l) for m in range(-l, l+1)]

            lmmax = len(LL)
            self.paw_qij = np.zeros((lmmax, lmmax))
            for ii in range(lmmax):
                for jj in range(ii+1):
                    if (LL[ii] == LL[jj]) and (MM[ii] == MM[jj]):
                        self.paw_qij[ii, jj] = self.radial_simp_int(
                            self.paw_ae_wfc[WW[ii]] * self.paw_ae_wfc[WW[jj]] -
                            self.paw_ps_wfc[WW[ii]] * self.paw_ps_wfc[WW[jj]]
                        )
                        self.paw_qij[jj, ii] = self.paw_qij[ii, jj]

        return self.paw_qij

    @property
    def symbol(self):
        '''
        return the symbol of the element
        '''
        return self.element

    @property
    def lmax(self):
        '''
        Return total number of l-channel projector functions.
        '''

        return self.proj_l.size

    @property
    def lmmax(self):
        '''
        Return total number of lm-channel projector functions.
        '''

        return np.sum(2 * self.proj_l + 1)

    def plot(self):
        '''
        '''
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        mpl.rcParams['axes.unicode_minus'] = False
        plt.style.use('ggplot')

        figure = plt.figure(
            figsize=(6.4, 6.4),
            # figsize = plt.figaspect(0.6),
            dpi=240,
        )

        axes = [
            plt.subplot(221),   # for real space projector functions
            plt.subplot(222),   # for reciprocal space projector functions
            plt.subplot(212)    # for ps/ae partial waves
        ]

        for ii in range(self.lmax):
            axes[0].plot(
                # self.proj_rgrid, self.rprojs[ii], label=f"L = {self.proj_l[ii]}"
                self.proj_rgrid, self.rprojs[ii], label="L = {}".format(
                    self.proj_l[ii])
            )
            axes[1].plot(
                # self.proj_qgrid, self.qprojs[ii], label=f"L = {self.proj_l[ii]}"
                self.proj_qgrid, self.qprojs[ii], label="L = {}".format(
                    self.proj_l[ii])
            )

            l1, = axes[2].plot(
                # self.rgrid, self.paw_ae_wfc[ii], label=f"L = {self.proj_l[ii]}"
                self.rgrid, self.paw_ae_wfc[ii], label=r"$L = {}$".format(
                    self.proj_l[ii])
            )
            axes[2].plot(
                self.rgrid, self.paw_ps_wfc[ii], ls=':',
                color=l1.get_color(),
                # label=r"$\tilde\phi_l(r)$ l = {}".format(self.proj_l[ii])
            )

        for ax in [axes[0], axes[2]]:
            ax.set_xlabel(r'$r\ [\AA]$', labelpad=5)

            ax.axhline(y=0, ls=':', color='k', alpha=0.6)
            ax.axvline(x=self.proj_rmax, ls=':', color='k', alpha=0.6)

        axes[2].legend(loc='best', fontsize='small', ncol=2)

        axes[1].set_xlabel(r'$G\ [\AA^{-1}]$', labelpad=5)
        axes[1].axhline(y=0, ls=':', color='k', alpha=0.6)
        axes[1].axvline(x=self.proj_gmax, ls=':', color='k', alpha=0.6)

        axes[0].set_title("Real-space Projectors", fontsize='small')
        axes[1].set_title("Reciprocal-space Projectors", fontsize='small')
        axes[2].set_title("AE/PS Partial Waves", fontsize='small')

        plt.tight_layout()
        plt.show()

    def __str__(self):
        '''
        '''
        # pstr = f"{self.symbol:>3s}\n"
        # pstr += f"\n{'l':>3s}{'rmax':>12s}\n"
        # pstr += ''.join(
        #     [f"{self.proj_l[ii]:>3d}{self.proj_rmax:>12.6f}\n"
        #         for ii in range(self.lmax)]
        # )

        pstr = "{:>3s}\n".format(self.symbol)
        pstr += "\n{:>3s}{:>12s}\n".format("l", "rmax")
        pstr += ''.join(
            ["{:>3d}{:>12.6f}\n".format(self.proj_l[ii], self.proj_rmax)
                for ii in range(self.lmax)]
        )
        return pstr


class nonlr(object):
    '''
    Real space presentation of the nonlocal projector functions on a regular 3d grid.

        p_{l,m}(r - R) = sqrt(Omega) * f(r - R) * ylm(r - R) * exp(ik*(r - R))

    where "f(r - R)" is the radial part of the real-space projector functions,
    which are stored in POTCAR file. "ylm(r - R)" is the real spherical harmonics
    with corresponding "l" and "m". "Omega" is the volume of the cell. The phase
    factor of "exp(ik*(r - R))" is stored in "crrexp".

    The application of the projector functions on the pseudo-wavefunction can
    then be obtained: C_n = < p_{l,m}(r - R) | \phi_{n,k} >

        C_n = \sum_{|r - R| < rmax} phi(r) * p_{l,m}(r - R)

    '''

    def __init__(self,
                 atoms, encut, potcar='POTCAR', k=[0.0, 0.0, 0.0],
                 ngrid=None,
                 lgam=False, lsoc=False
                 ):
        '''
        input:
            atoms: ase atom object
            encut: float, energy cutoff in eV
            potcar: the PAW POTCAR file of all the elements in atoms
            k: the k-point vector in fractional coordinate
        '''
        self.atoms = atoms
        self.natoms = len(atoms)
        self.encut = encut
        self.kvec = np.asarray(k, dtype=float)
        self.pawpp = [pawpotcar(potstr) for potstr in
                      open(potcar).read().split('End of Dataset')[:-1]]
        elements, elem_first_idx, elem_cnts = np.unique(atoms.get_chemical_symbols(),
                                                        return_index=True,
                                                        return_counts=True)
        # Sometimes, the order of the elements returned by np.unique may not be
        # consistent with that in POSCAR/POTCAR
        elem_first_idx = np.argsort(elem_first_idx)
        elements = elements[elem_first_idx]
        self.elem_cnts = elem_cnts[elem_first_idx]

        assert len(self.elem_cnts) == len(self.pawpp), \
            "The number of elements in POTCAR and POSCAR does not match!"
        if not np.alltrue([
            self.pawpp[ii].element.split('_')[0] == elements[ii]
            for ii in range(len(elements))
        ]):
            print(
                "\nWARNING:\nThe name of elements in POTCAR and POSCAR does not match!\n\n" +
                "    POTCAR: {}\n".format(' '.join([pp.element for pp in self.pawpp])) +
                "    POSCAR: {}\n".format(' '.join(elements))
            )

        self.elements = list(elements)
        self.element_idx = [self.elements.index(s) for s in
                            atoms.get_chemical_symbols()]

        self.set_grid(ngrid=ngrid)
        self.rphase()
        self.calc_rproj()

    def set_grid(self, ngrid):
        '''
        Minimum FFT grid size
        '''

        # real space cell
        self.Acell = self.atoms.get_cell()
        self.Anorm = np.linalg.norm(self.Acell, axis=1)
        # reciprocal space cell
        self.Bcell = np.linalg.inv(self.Acell).T
        self.Bnorm = np.linalg.norm(self.Bcell, axis=1)

        if ngrid is None:
            CUTOF = np.ceil(
                np.sqrt(self.encut / RYTOEV) / (TPI / (self.Anorm / AUTOA))
            )
            ################################################################################
            # In order to compare with VASP Normalcar, the grid size must be exactly
            # the same!
            ################################################################################
            self._ngrid = fftchk(np.array(2 * CUTOF, dtype=int) * 2)
        else:
            self._ngrid = np.array(ngrid, dtype=int)

    def rphase(self):
        '''
        Find lattice points contained within the PAW cutoff-sphere for each
        ion.
        '''
        self.ion_grid_idx = []
        self.ion_grid_direction = []
        self.ion_grid_distance = []
        self.ion_crrexp = []

        scaled_positions = self.atoms.get_scaled_positions()
        for iatom in range(self.natoms):
            # the type of the ion
            ntype = self.element_idx[iatom]
            # the paw pp of the ion
            pp = self.pawpp[ntype]
            # the positon of the ion in fractional coordinate
            R0 = scaled_positions[iatom]
            # PAW cutoff sphere radius of the ion
            rmax = pp.proj_rmax

            # restrict to points contained within a cube around the ion
            # grid position of the ion
            N0 = np.array(R0 * self._ngrid, dtype=int)
            # length of the cube
            ND = np.array(rmax * self.Bnorm * self._ngrid, dtype=int) + 1
            Dxyz = np.mgrid[N0[0] - ND[0]: N0[0] + ND[0] + 1,
                            N0[1] - ND[1]: N0[1] + ND[1] + 1,
                            N0[2] - ND[2]: N0[2] + ND[2] + 1].reshape((3, -1)).T
            # np.savetxt('dxyz.txt', Dxyz, fmt="%5d")
            # print(self._ngrid, Dxyz.max(),  Dxyz.min())

            # distance from the grid to the ion position: r - R0, fractional
            # coordinate
            # rR_d = Dxyz # / self._ngrid.astype(float)
            # print(rR_d.max(), rR_d.min())

            # distance from the grid to the ion position: r - R0
            rR = np.dot((Dxyz / self._ngrid.astype(float)) - R0, self.Acell)
            # print(rR.max(), rR.min(), rmax, self.Anorm)
            rRLen = np.linalg.norm(rR, axis=1)
            grid_in_sphere = rRLen <= (rmax / pp.NPSRNL * (pp.NPSRNL - 1))
            # print(rRLen[grid_in_sphere].max(), rmax)

            self.ion_grid_idx.append(Dxyz[grid_in_sphere] % self._ngrid)
            self.ion_grid_direction.append(rR[grid_in_sphere])
            # np.savetxt('rR{:d}'.format(iatom), rR[grid_in_sphere], fmt='%8.4f')
            self.ion_grid_distance.append(rRLen[grid_in_sphere])
            self.ion_crrexp.append(
                np.exp(1j * TPI * np.sum(
                    rR[grid_in_sphere] * np.dot(self.kvec, self.Bcell),
                    axis=1)))

    def calc_rproj(self):
        '''
        Nonlocal projector for each elements
        '''

        self.rproj_atoms = []
        for iatom in range(self.natoms):
            ntype = self.element_idx[iatom]
            pp = self.pawpp[ntype]
            # number of grid points in the sphere
            irmax = self.ion_grid_distance[iatom].shape[0]

            tmp = np.zeros((pp.lmmax, irmax))
            rproj_ylm = [
                sph_r(self.ion_grid_direction[iatom], l).T
                for l in range(pp.proj_l.max()+1)
            ]

            # np.savetxt('pj.{}'.format(pp.symbol), np.c_[pp.proj_rgrid, pp.rprojs.T])
            # xxx = np.zeros((pp.lmmax + 1, irmax))
            # xxx[0] = self.ion_grid_distance[iatom]
            # ii = 1

            iL = 0
            for l, spl_r in zip(pp.proj_l, pp.spl_rproj):
                TLP1 = 2 * l + 1
                rproj_radial = spl_r(self.ion_grid_distance[iatom])
                tmp[iL:iL+TLP1, :] = rproj_radial * rproj_ylm[l]
                iL += TLP1

            #     xxx[ii] = rproj_radial
            #     ii += 1
            # np.savetxt('pj.csp.{}{}'.format(pp.symbol, iatom), xxx.T)

            # For reciprocal space projectors, the factor is sqrt(1 / Omega)
            tmp *= np.sqrt(self.atoms.get_volume())
            self.rproj_atoms.append(tmp)

    def proj(self, wfc_r, whichatom=None):
        '''
        '''
        wfc_r = np.asarray(wfc_r)
        assert np.allclose(
            wfc_r.shape, self._ngrid), "Grid size does not match!"

        if whichatom is None:
            beta = []
            for iatom in range(self.natoms):
                gidx = self.ion_grid_idx[iatom]
                beta += [x for x in np.sum(
                    wfc_r[gidx[:, 0], gidx[:, 1], gidx[:, 2]] *
                    self.ion_crrexp[iatom] * self.rproj_atoms[iatom],
                    axis=1)
                ]
        else:
            gidx = self.ion_grid_idx[whichatom]
            beta = [x for x in np.sum(
                wfc_r[gidx[:, 0], gidx[:, 1], gidx[:, 2]] *
                self.ion_crrexp[whichatom] * self.rproj_atoms[whichatom])
            ]

        return np.asarray(beta)


class nonlq(object):
    '''
    Reciprocal space presentation of the nonlocal projector functions on a
    plane-waves grid.

        p_{l,m; R}(G + k) = (1. / sqrt(Omega) ) *\
                            i^l * f(G + k) * ylm(G + k) * exp(i(G+k)*R)

    where "f(G + k)" is the radial part of the reciprocal projector functions,
    which are stored in POTCAR file. "ylm(G+k)" is the real spherical harmonics
    with corresponding "l" and "m". "Omega" is the volume of the cell. The phase
    factor of "exp(i(G+k)*R)" is stored in "crexp". The "i^l" is stored in
    "cqfak".

    Note that the phase factor derives from the **translation/shifting
    property** of the FT.

    [Note 2021-11-28] Refer to this link for radial function FT.
    http://staff.ustc.edu.cn/~zqj/posts/FourierTransform-Radial-Function/

    The application of the projector functions on the pseudo-wavefunction can
    then be obtained: C_n = < p_{l,m; R} | \phi_{n,k} >

        C_n = \sum_G C_{n,k}(G + k) * p_{l,m}(G + k)

    '''

    def __init__(self,
                 atoms, encut, potcar='POTCAR', k=[0.0, 0.0, 0.0],
                 lgam=False, lsoc=False
                 ):
        '''
        input:
            atoms: ase atom object
            encut: float, energy cutoff in eV
            potcar: the PAW POTCAR file of all the elements in atoms
            k: the k-point vector in fractional coordinate
        '''
        self.atoms = atoms
        self.natoms = len(atoms)
        self.kvec = np.asarray(k, dtype=float)
        self.pawpp = [pawpotcar(potstr) for potstr in
                      open(potcar).read().split('End of Dataset')[:-1]]
        elements, elem_first_idx, elem_cnts = np.unique(atoms.get_chemical_symbols(),
                                                        return_index=True,
                                                        return_counts=True)
        # Sometimes, the order of the elements returned by np.unique may not be
        # consistent with that in POSCAR/POTCAR
        elem_first_idx = np.argsort(elem_first_idx)
        elements = elements[elem_first_idx]
        self.elem_cnts = elem_cnts[elem_first_idx]

        assert len(self.elem_cnts) == len(self.pawpp), \
            "The number of elements in POTCAR and POSCAR does not match!"
        if not np.alltrue([
            self.pawpp[ii].element.split('_')[0] == elements[ii]
            for ii in range(len(elements))
        ]):
            print(
                "\nWARNING:\nThe name of elements in POTCAR and POSCAR does not match!\n\n" +
                "    POTCAR: {}\n".format(' '.join([pp.element for pp in self.pawpp])) +
                "    POSCAR: {}\n".format(' '.join(elements))
            )

        self.elements = list(elements)
        self.element_idx = [self.elements.index(s) for s in
                            atoms.get_chemical_symbols()]
        # G-vectors in fractional coordinate
        self.Gvec = gvectors(atoms.cell, encut, k)
        self.nplw = self.Gvec.shape[0]
        # (k + G)-vectors in Cartesian coordinate
        self.Gk = np.dot(
            self.Gvec + self.kvec, TPI * self.atoms.get_reciprocal_cell()
        )
        # G-vectors length
        self.Glen = np.linalg.norm(self.Gk, axis=1)

        #
        self.setylm()
        self.phase()
        self.calc_qproj()

    def setylm(self):
        '''
         Calculate the real spherical harmonics for a set of G-grid points up to
         LMAX.
        '''

        lmax = np.max([p.proj_l.max() for p in self.pawpp])
        self.ylm = []
        for l in range(lmax+1):
            self.ylm.append(
                sph_r(self.Gk, l)
            )

    def phase(self):
        '''
        Calculates the phasefactor CREXP (exp(i(G + k).R)) for one k-point
        '''
        #####################################################################
        # Mind the sigh of "1j" here. I used "-1j" at first, which took me a
        # long long time to figure out what goes wrong!
        #####################################################################
        self.crexp = np.exp(1j * TPI *
                            np.dot(
                                self.Gvec, self.atoms.get_scaled_positions().T
                            ))
        # i^L is stored in CQFAK
        self.cqfak = [
            1j ** np.array([
                l for l in pp.proj_l
                for ii in range(2 * l + 1)
            ])
            for pp in self.pawpp
        ]

    def calc_qproj(self):
        '''
        Nonlocal projector for each type of element.
        '''
        self.qproj = []
        for pp in self.pawpp:
            # np.savetxt('pj.{}'.format(pp.symbol), np.c_[pp.proj_qgrid, pp.qprojs.T])
            # xxx = np.zeros((pp.lmmax + 1, self.nplw))
            # xxx[0] = self.Glen
            # ii = 1

            # find out those | G + k | <= gmax of reciprocal projectors
            G_within_gmax = (self.Glen <= pp.proj_gmax)

            tmp = np.zeros((pp.lmmax, self.nplw))
            iL = 0
            for l, spl_q in zip(pp.proj_l, pp.spl_qproj):
                TLP1 = 2 * l + 1
                # radial part of the projector: spl_q(self.Glen)
                # spherical harmonics of angular momentum l: ylm[l]
                tmp[iL:iL+TLP1, G_within_gmax] = spl_q(self.Glen[G_within_gmax]) *\
                    self.ylm[l].T[:, G_within_gmax]

                iL += TLP1

                # xxx[ii] = spl_q(self.Glen)
                # ii += 1
            # np.savetxt('pj.csp.{}'.format(pp.symbol), xxx.T)

            tmp /= np.sqrt(self.atoms.get_volume())
            self.qproj.append(tmp)

    def proj(self, cptwf, whichatom=None):
        '''
        Project one single KS wavefunctions onto all the nonlocal reciprocal
        space projectors.
        '''

        cptwf = np.asarray(cptwf)
        assert cptwf.size == self.nplw, "Number of plane waves does not match!"

        if whichatom is None:
            beta = []
            for iatom in range(self.natoms):
                ntype = self.element_idx[iatom]
                iill = self.cqfak[ntype]
                beta += [x for x in
                         np.sum(
                             cptwf * self.crexp[:, iatom] *
                             (self.qproj[ntype] * iill[:, None]),
                             axis=1
                         )]
        else:
            ntype = self.element_idx[whichatom]
            iill = self.cqfak[ntype]
            beta = [x for x in
                    np.sum(
                        cptwf * self.crexp[:, whichatom] *
                        (self.qproj[ntype] * iill[:, None]),
                        axis=1
                    )]
        return np.asarray(beta)


class radial2grid(object):
    '''
    '''

    def __init__(self,
                 r, fr, cell, encut,
                 R0=[0.0, 0.0, 0.0],
                 # bc_type='natural',
                 rlog=False,
                 reciprocal=False):
        '''
        inputs
            r: the coordinate of the radial grid (r-grid)
            fr: the function values on the r-grid
            cell:  (3,3) ndarray in units of Angstrom, the basis vectors of the regular grid
            encut: the energy cutoff, which determines the grid size of the
                   regurlar grid
            R0: coordinate of the center of the core region
            # bc_type: boundary condition to interpolate r/fr on the radial grid
            rlog: logarithmic radial grid?
            reciprocal: r/fr defined in reciprocal space?
        '''
        from scipy.interpolate import CubicSpline as csp

        if not rlog:
            if reciprocal:
                fr_cs = csp(r, fr, bc_type='natural')
            else:
                fr_cs = csp(r, fr, bc_type='natural')
        else:
            pass


if __name__ == '__main__':
    import time
    xx = open('examples/projectors/lreal_true/potcar.mo').read()

    t0 = time.time()
    ps = pawpotcar(xx)

    # t1 = time.time()
    # ps.csplines()
    # t2 = time.time()
    # print(t1 - t0)
    # print(t2 - t1)

    # print(ps.symbol)
    # print(ps.lmmax, ps.lmax)
    # print(ps)
    # print(ps.paw_ae_wfc[1][-1])

    ps.plot()
