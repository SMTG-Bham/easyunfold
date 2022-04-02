"""
Procar reader
"""
from typing import List, Union
from pathlib import Path
import re
import numpy as np


class Procar:
    """Reader for PROCAR file"""

    def __init__(self, fobj_or_path=None, is_soc=False):
        """
        Read the PROCAR file
        """
        self._is_soc = is_soc
        self.eigenvalues = None
        self.kvecs = None
        self.kweights = None
        self.nbands = None
        self.nkpts = None
        self.nspins = None
        self.nion = None
        self.occs = None
        self.proj_names = None
        self.proj_data = None
        self.header = None
        self.proj_xyz = None

        # Read the PROCAR
        if isinstance(fobj_or_path, (str, Path)):
            with open(fobj_or_path) as fhandle:
                self._read(fhandle)
        else:
            self._read(fobj_or_path)

    def _read(self, fobj):
        """Main function for reading in the data"""

        # First sweep - found the number of kpoints and the number of bands
        fobj.seek(0)
        self.header = fobj.readline()
        # Read the NK, NB and NIONS that are integers
        self.nkpts, self.nbands, self.nion = [int(token) for token in re.sub(r'[^0-9]', ' ', fobj.readline()).split()]
        # Number of projects and their names
        nproj = None
        self.proj_names = None
        for line in fobj:
            if re.match(r'^ion', line):
                nproj = len(line.strip().split()) - 2
                self.proj_names = line.strip().split()[1:-1]
                break
        # Count the number of data lines, these lines do not have any alphabets
        proj_data = []
        energies = []
        occs = []
        kvecs = []
        kweights = []
        fobj.seek(0)
        for line in fobj:
            if not re.search(r'[a-zA-Z]', line) and line.strip():
                proj_data.append([float(token) for token in line.strip().split()[1:-1]])
            elif line.startswith('band'):
                tokens = line.strip().split()
                energies.append(float(tokens[4]))
                occs.append(float(tokens[-1]))
            elif line.startswith(' k-point'):
                line = re.sub(r'(\d)-', r'\1 -', line)
                tokens = line.strip().split()
                kvecs.append([float(val) for val in tokens[-6:-3]])
                kweights.append(float(tokens[-1]))
        self.occs = np.array(occs)
        self.kvecs = np.array(kvecs)
        self.kweights = np.array(kweights)
        self.eigenvalues = np.array(energies)

        proj_data = np.array(proj_data, dtype=float)

        self.nspins = proj_data.shape[0] // (self.nion * self.nbands * self.nkpts)
        self.nspins //= 4 if self._is_soc else 1
        if self._is_soc:
            self.spin = self.spin // 4
        # Reshape
        self.occs.resize((self.nspins, self.nkpts, self.nbands))
        self.kvecs.resize((self.nspins, self.nkpts, 3))
        self.kweights.resize((self.nspins, self.nkpts))
        self.eigenvalues.resize((self.nspins, self.nkpts, self.nbands))

        # Reshape the array
        if self._is_soc is False:
            self.proj_data = proj_data.reshape((self.nspins, self.nkpts, self.nbands, self.nion, nproj))
        else:
            self.proj_data = proj_data.reshape((self.nspins, self.nkpts, self.nbands, 4, self.nion, nproj))
            # Split the data into xyz projection and total
            self.proj_xyz = self.proj_data[:, :, :, 1:, :, :]
            self.proj_data = self.proj_data[:, :, :, 0, :, :]

    def get_projection(self, atom_idx: List[int], proj: Union[List[str], str], weight_by_k=False):
        """
        Get project for specific atoms and specific projectors

        Args:
            atom_idx (list): A list of index of the atoms to be selected
            proj (list): A list of the projector names to be selected
            weight_by_k: Apply k weighting or not.
        Returns:
            The project summed over the selected atoms and the projectors
        """
        atom_mask = [iatom in atom_idx for iatom in range(self.nion)]
        assert any(atom_mask)
        if proj == 'all':
            out = self.proj_data[:, :, :, atom_mask, :].sum(axis=(-1, -2))
        else:
            proj_idx = [proj_name in proj for proj_name in self.proj_names]
            assert any(proj_idx)
            out = self.proj_data[:, :, :, :, proj_idx]
            out = out[:, :, :, atom_idx, :].sum(axis=(-1, -2))

        if weight_by_k:
            for kidx in range(self.nkpts):
                out[:, kidx, :] *= self.kweights[kidx]
        return out
