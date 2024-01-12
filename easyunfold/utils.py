"""
Utility functions
"""
import re
from pathlib import Path
from typing import Union, List
import numpy as np
from castepinput import CellInput, Block

RE_COMMENT = re.compile(r'[!#]')


def write_kpoints(kpoints: Union[np.ndarray, list], outpath, *args, code='vasp', **kwargs):
    """
    Write the kpoints to a file

    :param kpoints: A Nx3 array of the kpoint coordinates (fractional) to be written.
    :param outpath: Path of the output file to be used.
    :param code: The code that the kpoint file is used for.
    """
    kpoints = np.asarray(kpoints)

    if code == 'vasp':
        return write_kpoints_vasp(kpoints, outpath, *args, **kwargs)
    if code == 'castep':
        return write_kpoints_castep(kpoints, outpath, *args, **kwargs)
    raise NotImplementedError(f'DFT code: {code} is not implemented')


def write_kpoints_castep(kpoints: Union[np.ndarray, list], dest, source=None, tag='spectral', weights=None, **kwargs):
    """
    Write kpoints to a CASTEP input file.
    Optically can use an existing file as the template and update it.
    """
    _ = kwargs
    if source is not None:
        cell = CellInput.from_file(source)
    else:
        cell = CellInput()

    if weights is None:
        cell[f'{tag}_kpoints_list'.upper()] = Block([f'{k[0]:.10f} {k[1]:.10f} {k[2]:.10f}' for k in kpoints])
    else:
        cell[f'{tag}_kpoints_list'.upper()] = Block([f'{k[0]:.10f} {k[1]:.10f} {k[2]:.10f} {w:.10f}' for k, w in zip(kpoints, weights)])

    cell.save(dest)


def write_kpoints_vasp(kpoints: Union[np.ndarray, list],
                       outpath: str = 'KPOINTS',
                       comment: str = '',
                       weights: Union[None, List[float]] = None,
                       **kwargs):
    """
    Write kpoints to VASP KPOINTS file

    :param kpoints: A list of kpoints to be written
    :param outpath: Path to the output file
    :param comments: Comments to be put into the `KPOINTS` file
    :param weights: If given, the weighting of the kpoints, otherwise all kpoints are equal weighted.
    """
    _ = kwargs
    kpoints = np.asarray(kpoints)
    nkpts = kpoints.shape[0]

    with open(outpath, 'w', encoding='utf-8') as vaspkpt:
        vaspkpt.write(comment + '\n')
        vaspkpt.write(f'{nkpts}\n')
        vaspkpt.write('Rec\n')
        for ik in range(nkpts):
            if weights is not None:
                line = f'  {kpoints[ik, 0]:12.8f} {kpoints[ik, 1]:12.8f} {kpoints[ik,2]:12.8f}  {weights[ik]:12.8f}\n'
            else:
                line = f'  {kpoints[ik, 0]:12.8f} {kpoints[ik, 1]:12.8f} {kpoints[ik,2]:12.8f}  1.0\n'
            vaspkpt.write(line)


def read_kpoints_line_vasp(content, density=20):
    """
    Read kpoints in the VASP KPOINTS file line mode, the results are resolved
    to explicit kpoints

    :returns: The kpoints, the comment, the labels, and the weights at each kpoint
    """
    comment = content[0]
    density = int(content[1]) if content[1] else density

    assert content[2].lower().startswith('lin'), 'Only Line mode KPOINT file is supported'
    assert content[3].lower().startswith('rec'), 'Only Reciprocal coorindates are supported!'

    segs = []
    labels = []
    for line in content[4:]:
        tokens = line.split()
        if not tokens:
            continue
        point = [float(x) for x in tokens[:3]]
        labels.append(tokens[-1])
        segs.append(point)
    # Process the segments
    kpoints = []
    labels_loc = []
    for i in range(int(len(segs) / 2)):
        k1 = segs[i * 2]
        k2 = segs[i * 2 + 1]
        # Check for duplicate end point
        if kpoints and kpoints[-1] == k1:
            kpoints.pop()
            labels_loc.pop()
        this_seg = np.linspace(k1, k2, density)
        labels_loc.append((len(kpoints), labels[i]))
        kpoints.extend(this_seg.tolist())
        labels_loc.append((len(kpoints) - 1, labels[i + 1]))
    return kpoints, comment, labels_loc, None


def read_kpoints(path='KPOINTS', code='vasp', **kwargs):
    """
    Read the kpoints from a given file

    This function dispatches to code-specific function based on the `code` keyword variable.

    """
    if code == 'vasp':
        return read_kpoints_vasp(path, **kwargs)
    if code == 'castep':
        return read_kpoints_castep(path, **kwargs)
    raise NotImplementedError(f'DFT code: {code} is not implemented')


def read_kpoints_vasp(path='KPOINTS'):
    """
    Read kpoints from a KPOINTS file containing reciprocal space coordinates (fractional)

    :returns: The kpoints, the comment, the labels, and the weights at each kpoint.
    """
    content = Path(path).read_text(encoding='utf-8').split('\n')
    comment = content[0]
    nkpts = int(content[1])
    if content[2].lower().startswith('lin'):
        return read_kpoints_line_vasp(content)
    assert content[2].lower().startswith('rec'), 'Only Reciprocal space KPOINT file is supported'
    kpts = []
    labels = []
    weights = []
    ik = 0
    for line in content[3:]:
        tokens = line.split()
        this_kpt = [float(value) for value in tokens[:3]]
        weights.append(float(tokens[3]))

        if len(tokens) >= 5:
            labels.append([ik, tokens[4]])

        kpts.append(this_kpt)
        ik += 1
        if ik == nkpts:
            break
    return kpts, comment, labels, weights


def read_kpoints_castep(path, tag='spectral'):
    """
    Read explicitly defined kpoints in a cell file to CASTEP
    """
    lines = Path(path).read_text(encoding='utf-8').split('\n')
    tags = [f'{tag}_kpoints_list', f'{tag}_kpoint_list']
    capture = False
    kpts = []
    weights = []
    labels = []
    for line in lines:
        # Skip and empty lines
        if not line.strip() or line.startswith('#'):
            continue

        if r'%block' in line.lower():
            block_name = line.split()[1].lower()
            if block_name in tags:
                capture = True
                continue

        if r'%endblock' in line.lower() and capture:
            block_name = line.split()[1].lower()
            if block_name not in tags:
                raise RuntimeError(f'Unexpected end of block detected {line}')
            break

        # Capturing mode - we are inside of a valid block
        if capture:
            tokens = line.strip().split()
            kpts.append(list(float(tmp) for tmp in tokens[:3]))
            if len(tokens) > 3:
                # Check if the weights are included
                if not RE_COMMENT.match(tokens[3]):
                    weights.append(float(tokens[3]))
            sub_lines = RE_COMMENT.split(line)
            if len(sub_lines) > 1:
                # Handle case like: 0 0 0 0.1 # \Gamma
                # Record the label and the index of the kpoint
                labels.append([len(kpts) - 1, sub_lines[1].split()[0]])
    if not weights:
        weights = None
    return kpts, '', labels, weights


def wrap_kpoints(kpoints: Union[list, np.ndarray]):
    """Wrap the kpoints to range [-0.5, 0.5)"""
    kpoints = np.array(kpoints) + 0.5
    kpoints -= np.floor(kpoints)
    kpoints -= 0.5
    # Giving some numerical tolerance when enforcing the range [-0.5, 0.5)
    kpoints[np.abs(kpoints - 0.5) < 1e-7] = -0.5
    return kpoints


def find_unique(seq: np.ndarray, func=None):
    """
    Find unique slices along the first dimension of an np.array.
    This function is not optimised for high performance and has a O(N^2) scaling.

    :returns: A tuple of (unique, unique_idx, inv_mapping)
    """
    if func is None:
        # Use equality condition
        def _func(x, y):
            """Check elements of x and y are all the same"""
            return np.all(x == y)

        func = _func

    mapping_inv = np.zeros(len(seq), dtype=int) - 1
    unique_idx = []
    nobj = len(seq)
    for i in range(nobj):
        # Have this object been mapped?
        if mapping_inv[i] >= 0:
            continue
        # This object has not been mapped to any unique obj identified
        unique_idx.append(i)
        mapping_inv[i] = len(unique_idx) - 1
        # Forward search for any object that is identical with this object
        for j in range(i + 1, nobj):
            if func(seq[i], seq[j]):
                # j is the same as i
                mapping_inv[j] = len(unique_idx) - 1

    unique_idx = np.array(unique_idx)
    return seq[unique_idx], unique_idx, mapping_inv


def reduce_kpoints(kpoints: Union[list, np.ndarray], time_reversal=True, rounding_digits=10):
    """
    Reduce the kpoint set by finding duplicated kpoints
    """
    kpoints = np.asarray(kpoints)
    kpoints_rounded = np.round(wrap_kpoints(kpoints), rounding_digits)

    if not time_reversal:
        # No time-reversal - use np.unique for speed
        _, unique_id, inv_mapping = np.unique(kpoints_rounded, axis=0, return_inverse=True, return_index=True)
    else:

        def equality_time_reversal(x, y):
            """Check if x == y or x == -y"""
            return np.all(x == y) | np.all(x == -y)

        _, unique_id, inv_mapping = find_unique(kpoints_rounded, equality_time_reversal)

    unique_k = kpoints[unique_id]
    return unique_k, unique_id, inv_mapping


def kpoints_equal(k1, k2, time_reversal=False, atol=1e-5):
    """Return two if two kpoints are equivalent to each other"""
    if np.allclose(wrap_kpoints(k1), wrap_kpoints(k2), atol=atol):
        return True
    if time_reversal and np.allclose(wrap_kpoints(k1), -wrap_kpoints(k2), atol=atol):
        return True
    return False
