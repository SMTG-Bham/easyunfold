"""
Procar reader
"""

import os.path
from typing import List, Union
from pathlib import Path
import re
import numpy as np

from monty.io import zopen
from monty.json import MSONable, MontyDecoder

from easyunfold import __version__

# pylint:disable=too-many-locals,


class Procar(MSONable):
    """Reader for PROCAR file"""

    def __init__(self, fobjs_or_paths=None, is_soc=False, normalise=True):
        """
        Read the PROCAR(.gz) file from a handle or path

        :param fobjs_or_paths:  Either a string or list of file-like objs or paths
        :param is_soc: Whether the PROCAR(.gz) is from a calculation with spin-orbit coupling
        :param normalise: Whether to normalise the projection for every band or not
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
        self.normalise = normalise

        # Read the PROCAR
        if isinstance(fobjs_or_paths, (str, Path)):
            fobjs_or_paths = [fobjs_or_paths]
        self.read(fobjs_or_paths)

    def _read(self, fobj, parsed_kpoints=None):
        """Main function for reading in the data"""
        if parsed_kpoints is None:
            parsed_kpoints = set()

        # First sweep - find the number of kpoints and the number of bands
        fobj.seek(0)
        _header = fobj.readline()
        # Read the NK, NB and NIONS that are integers
        _total_nkpts, nbands, nion = [int(token) for token in re.sub(r'[^0-9]', ' ', fobj.readline()).split()]
        if nion != self.nion:
            raise ValueError(f'Mismatch in number of ions in PROCARs supplied: ({nion} vs {self.nion})!')

        # Count the number of data lines, these lines do not have any alphabets
        proj_data, energies, kvecs, kweights, occs = [], [], [], [], []
        tot_count = 0  # count the instances of lines starting with "tot" -> (4 + 1) * nbands * nkpts for SOC calcs
        fobj.seek(0)

        line = fobj.readline()
        # Counter for the number of sections in the PROCAR
        section_counter = 1
        _last_kid = 0
        this_procar_parsed_kpoints = set()  # set with tuples of parsed (kvec tuple, section_counter) for this PROCAR
        while line:
            if line.startswith(' k-point'):
                line = re.sub(r'(\d)-', r'\1 -', line)
                tokens = line.strip().split()
                _kid = int(tokens[1])

                # Check if the VASP PROCAR k index decreases – if so we have entered a second section
                if _kid < _last_kid:
                    section_counter += 1
                _last_kid = _kid

                kvec = tuple(round(float(val), 5) for val in tokens[-6:-3]  # tuple to make it hashable
                            )  # round to 5 decimal places to ensure proper kpoint matching
                if (kvec not in parsed_kpoints and (kvec, section_counter) not in this_procar_parsed_kpoints):
                    this_procar_parsed_kpoints.add((kvec, section_counter))
                    kvecs.append(list(kvec))
                    kweights.append(float(tokens[-1]))
                else:
                    # skip ahead to the next instance of two blank lines in a row
                    while line.strip() or fobj.readline().strip():
                        line = fobj.readline()
                    continue

            elif (not re.search(r'[a-zA-Z]', line) and line.strip() and len(line.strip().split()) - 2 == len(self.proj_names)):
                # only parse data if line is expected length, in case of LORBIT >= 12
                proj_data.append([float(token) for token in line.strip().split()[1:-1]])

            elif line.startswith('band'):
                tokens = line.strip().split()
                energies.append(float(tokens[4]))
                occs.append(float(tokens[-1]))

            elif line.startswith('tot'):
                tot_count += 1

            line = fobj.readline()

        # dynamically determine whether PROCARs are SOC or not
        if tot_count == 4 * len(occs):
            self._is_soc = True
        elif tot_count == len(occs):
            self._is_soc = False
        else:
            raise ValueError(f"Number of lines starting with 'tot' ({tot_count}) in PROCAR does not match expected "
                             f'values ({4*len(occs)} or {len(occs)})!')

        occs = np.array(occs)
        kvecs = np.array(kvecs)
        kweights = np.array(kweights)
        eigenvalues = np.array(energies)

        proj_data = np.array(proj_data, dtype=float)

        nkpts = len(kvecs)  # redetermine nkpts in case some were skipped due to already being parsed

        # For spin-polarised calcs, the data from the second (down) spin are located after that of the first (up) spin
        # Hence, the number of spins is simply the number of sections
        self.nspins = section_counter

        # Reshape
        occs.resize((self.nspins, nkpts // self.nspins, nbands))
        kvecs.resize((self.nspins, nkpts // self.nspins, 3))
        kweights.resize((self.nspins, nkpts // self.nspins))
        eigenvalues.resize((self.nspins, nkpts // self.nspins, nbands))

        # Reshape the array
        if self._is_soc is False:
            proj_data = proj_data.reshape((
                self.nspins,
                nkpts // self.nspins,
                nbands,
                self.nion,
                len(self.proj_names),
            ))
            proj_xyz = None
        else:
            proj_data = proj_data.reshape((self.nspins, nkpts, nbands, 4, self.nion, len(self.proj_names)))
            # Split the data into xyz projection and total
            proj_xyz = proj_data[:, :, :, 1:, :, :]
            proj_data = proj_data[:, :, :, 0, :, :]

        if self.normalise:
            self.normalise_projs(proj_data)

            if proj_xyz is not None:
                proj_sum = np.sum(proj_xyz, axis=(-3, -2, -1), keepdims=True)
                proj_sum[proj_sum == 0] = 1
                proj_xyz /= proj_sum

        # Update the parsed kpoints
        parsed_kpoints.update({kvec_section_counter_tuple[0] for kvec_section_counter_tuple in this_procar_parsed_kpoints})

        return (
            self.nspins,
            occs,
            kvecs,
            kweights,
            eigenvalues,
            proj_data,
            proj_xyz,
            parsed_kpoints,
        )

    def normalise_projs(self, proj_data):
        """
        Normalise the projections

        For each nspin, nkpt, nband, normalise the sum of projections over nion and proj_names to be 1.
        Atomic & orbital projections do not sum to 1 in most cases in VASP, as only those falling inside
        the atomic radii and overlapping with spd spherical harmonics are counted.
        """
        proj_sum = np.sum(proj_data, axis=(-2, -1), keepdims=True)
        proj_sum[proj_sum == 0] = 1  # just in case, avoid division by zero
        proj_data /= proj_sum
        self.normalise = True

    def _read_header_nion_proj_names(self, fobj):
        """Read the header, nion and proj_names from the PROCAR"""
        fobj.seek(0)
        self.header = fobj.readline()
        # Read the NK, NB and NIONS that are integers
        _nkpts, _nbands, self.nion = [int(token) for token in re.sub(r'[^0-9]', ' ', fobj.readline()).split()]
        self.proj_names = None  # projection names

        for line in fobj:
            if re.match(r'^ion.*tot', line):  # only the first "ion" line, in case of LORBIT >= 12
                self.proj_names = line.strip().split()[1:-1]
                break

    def read(self, fobjs_or_paths):
        """Read and amalgamate the data from a list of PROCARs"""

        def open_file(fobj_or_path):
            if isinstance(fobj_or_path, (str, Path)):
                if os.path.exists(fobj_or_path):
                    return zopen(fobj_or_path, mode='rt')  # closed later
                if os.path.exists(f'{fobj_or_path}.gz'):
                    return zopen(f'{fobj_or_path}.gz', mode='rt')

                raise FileNotFoundError(  # else raise error
                    f'File not found: {fobj_or_path} – PROCAR(.gz) file needed for '
                    f'parsing atomic projections!')
            return fobj_or_path  # already a file-like object, just return it

        parsed_kpoints = None
        occs_list, kvecs_list, kweights_list = [], [], []
        eigenvalues_list, proj_data_list, proj_xyz_list = [], [], []
        for i, fobj_or_path in enumerate(fobjs_or_paths):
            # Note: If PROCAR parsing becomes a significant bottleneck for people (e.g. with several HSE06+SOC PROCARs),
            # this could be parallelized (somewhat) with multiprocessing. The actual file parsing in _read() is currently
            # serial so no easy wins there, but could at least parallelise over the list of PROCARs
            fobj = open_file(fobj_or_path)
            if self.header is None:  # first file; read header, nion, proj_names
                self._read_header_nion_proj_names(fobj)

            current_nspins = self.nspins  # check spin consistency between PROCARs
            (
                nspins,
                occs,
                kvecs,
                kweights,
                eigenvalues,
                proj_data,
                proj_xyz,
                parsed_kpoints,
            ) = self._read(fobj, parsed_kpoints=parsed_kpoints)
            if current_nspins is not None and current_nspins != nspins:
                raise ValueError(f'Mismatch in number of spins in PROCARs supplied: ({nspins} vs {current_nspins})!')

            if isinstance(fobj_or_path, (str, Path)):
                fobj.close()  # if file was opened in this loop, close it

            # Append to respective lists
            occs_list.append(occs)
            kvecs_list.append(kvecs)
            kweights_list.append(kweights)
            eigenvalues_list.append(eigenvalues)
            proj_data_list.append(proj_data)
            proj_xyz_list.append(proj_xyz)
            if len(fobjs_or_paths) > 1:  # print progress if reading multiple files
                print(f'Finished parsing PROCAR {i + 1}/{len(fobjs_or_paths)}')

        # Combine along the nkpts axis:
        # for occs, eigenvalues, proj_data and proj_xyz, nbands (axis = 2) could differ, so set missing values to zero:
        max_nbands = max(arr.shape[2] for arr in eigenvalues_list)
        for array_list in [occs_list, eigenvalues_list, proj_data_list, proj_xyz_list]:
            for i, arr in enumerate(array_list):
                if arr is not None and arr.shape[2] < max_nbands:
                    if len(arr.shape) == 3:  # occs_list, eigenvalues_list
                        array_list[i] = np.pad(
                            arr,
                            ((0, 0), (0, 0), (0, max_nbands - arr.shape[2])),
                            mode='constant',
                        )
                    elif len(arr.shape) == 5:  # proj_xyz_list
                        array_list[i] = np.pad(
                            arr,
                            (
                                (0, 0),
                                (0, 0),
                                (0, max_nbands - arr.shape[2]),
                                (0, 0),
                                (0, 0),
                            ),
                            mode='constant',
                        )
                    elif len(arr.shape) == 6:  # proj_xyz_list
                        array_list[i] = np.pad(
                            arr,
                            (
                                (0, 0),
                                (0, 0),
                                (0, max_nbands - arr.shape[2]),
                                (0, 0),
                                (0, 0),
                                (0, 0),
                            ),
                            mode='constant',
                        )
                    else:
                        raise ValueError('Unexpected array shape encountered!')

        self.nbands = max_nbands
        self.occs = np.concatenate(occs_list, axis=1)
        self.eigenvalues = np.concatenate(eigenvalues_list, axis=1)
        self.kvecs = np.concatenate(kvecs_list, axis=1)
        self.kweights = np.concatenate(kweights_list, axis=1)
        self.proj_data = np.concatenate(proj_data_list, axis=1)
        if all(arr is not None for arr in proj_xyz_list):
            self.proj_xyz = np.concatenate(proj_xyz_list, axis=1)

        self.nkpts = self.kvecs.shape[1]

    def get_projection(self, atom_idx: List[int], proj: Union[List[str], str], weight_by_k=False):
        """
        Get projection for specific atoms and specific projectors


        :param atom_idx: A list of index of the atoms to be selected
        :param proj: A list of the projector names to be selected
        :param weight_by_k: Apply k weighting or not.

        :returns: The projection summed over the selected atoms and the projectors
        """
        atom_mask = [iatom in atom_idx for iatom in range(self.nion)]
        assert any(atom_mask)
        if proj == 'all':
            out = self.proj_data[:, :, :, atom_mask, :].sum(axis=(-1, -2))
        else:
            if isinstance(proj, str):
                proj = [
                    proj,
                ]

            # replace any instance of "p" with "px,py,pz" and "d" with "dxy,dyz,dz2,dxz,dx2-y2"
            def _replace_p_d(single_proj):
                if single_proj == 'p':
                    return ['px', 'py', 'pz']
                if single_proj == 'd':
                    return [
                        'dxy',
                        'dyz',
                        'dz2',
                        'dxz',
                        'x2-y2',
                    ]  # dx2-y2 labelled differently in VASP
                    # PROCAR

                return [single_proj]

            proj = [_replace_p_d(single_proj) for single_proj in proj]
            proj = [item for sublist in proj for item in sublist]

            proj_idx = [proj_name in proj for proj_name in self.proj_names]
            assert any(proj_idx)
            out = self.proj_data[:, :, :, :, proj_idx]
            out = out[:, :, :, atom_idx, :].sum(axis=(-1, -2))

        if weight_by_k:
            for kidx in range(self.nkpts):
                out[:, kidx, :] *= self.kweights[:, kidx]
        return out

    def as_dict(self) -> dict:
        """Convert the object into a dictionary representation (so it can be saved to json)"""
        output = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': __version__,
        }
        for key in [
                '_is_soc',
                'eigenvalues',
                'kvecs',
                'kweights',
                'nbands',
                'nkpts',
                'nspins',
                'nion',
                'occs',
                'proj_names',
                'proj_data',
                'header',
                'proj_xyz',
                'normalise',
        ]:
            output[key] = getattr(self, key)
        return output

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs Procar object from a dict representation, without calling __init__().

        Args:
            d (dict): dict representation of Procar

        Returns:
            Procar object
        """

        def decode_dict(subdict):
            if isinstance(subdict, dict) and '@module' in subdict:
                return MontyDecoder().process_decoded(subdict)
            return subdict

        instance = cls.__new__(cls)  # create a new instance without calling __init__()
        d_decoded = {k: decode_dict(v) for k, v in d.items()}

        # set the instance variables directly from the dictionary
        for key, value in d_decoded.items():
            if key in ['@module', '@class', '@version']:
                continue
            setattr(instance, key, value)

        return instance
