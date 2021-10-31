# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import glob
import os
import struct
import numpy as np
from ..core import Array
from .. import config
from .. import units


def generate_fname(nout, path="", ftype="", cpuid=1, ext=""):

    if nout == -1:
        filelist = sorted(glob.glob(os.path.join(path, "output*")))
        number = filelist[-1].split("_")[-1]
    else:
        number = str(nout).zfill(5)

    infile = os.path.join(path, "output_" + number)
    if len(ftype) > 0:
        infile = os.path.join(infile, ftype + "_" + number)
        if cpuid > 0:
            infile += ".out" + str(cpuid).zfill(5)

    if len(ext) > 0:
        infile += ext

    return infile


def read_parameter_file(fname=None, delimiter="="):
    """
    Read info file and create dictionary
    """
    out = {}
    with open(fname, 'r') as f:
        content = f.readlines()
    for line in content:
        sp = line.split(delimiter)
        if len(sp) > 1:
            value = sp[1].strip()
            try:
                value = eval(value)
            except NameError:
                pass
            out[sp[0].strip()] = value
    return out


def read_binary_data(content=None,
                     fmt=None,
                     offsets=None,
                     skip_head=True,
                     increment=True):
    """
    Unpack binary data from a content buffer using a dict of offsets.
    Also increment the offsets of the corresponding data read, as well as
    increase the line count by 1.
    """

    byte_size = {
        "b": 1,
        "h": 2,
        "i": 4,
        "q": 8,
        "f": 4,
        "d": 8,
        "e": 8,
        "n": 8,
        "l": 8,
        "s": 1
    }

    offset = 0
    for key in offsets:
        offset += offsets[key] * byte_size[key]
    # if offset is None:
    #     offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16
    if skip_head:
        offset += 4  # + correction

    # byte_size = {"b": 1 , "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = int(fmt[:-1])
    pack_size = mult * byte_size[fmt[-1]]

    if increment:
        offsets[fmt[-1]] += mult
    offsets["n"] += 1

    return struct.unpack(fmt, content[offset:offset + pack_size])


def skip_binary_line(content, offsets):
    """
    Return the number of bytes necessary to skip the current line.
    """
    [nbytes] = read_binary_data(fmt="i",
                                content=content,
                                offsets=offsets,
                                skip_head=False,
                                increment=False)
    return nbytes


def make_vector_arrays(data, ndim):
    """
    Merge vector components in 2d arrays.
    """
    components = list("xyz"[:ndim])
    if len(components) > 1:
        delete = []
        for key in list(data.keys()):
            comp_list = None
            rawkey = None
            inds = [i for i, letter in enumerate(key) if letter == 'x']
            for ind in inds:
                comp_list = [key[:ind] + c + key[ind + 1:] for c in components]
                if all([item in data for item in comp_list]):
                    cut = ind - 1 if key[ind - 1] == "_" else ind
                    rawkey = key[:cut] + key[ind + 1:]
                    if len(rawkey) == 0:
                        rawkey = "xyz"
                    data[rawkey] = Array(values=np.array(
                        [data[c].values for c in comp_list]).T,
                                         unit=data[key].unit)
                    delete += comp_list
        for key in delete:
            del data[key]


def find_max_amr_level(levelmax, select):
    """
    Test the selection function in `select` on the range of possible AMR levels
    to determine the max level to read.
    """
    possible_levels = np.arange(1, levelmax + 1, dtype=int)
    func_test = select["level"](possible_levels)
    inds = np.argwhere(func_test).ravel()
    return possible_levels[inds.max()]


def get_spatial_scaling(ud, ul, ut, scale):
    """
    Compute the scaling factor to convert between code units and requested spatial
    scale.
    """
    length_unit = config.get_unit("x", ud, ul, ut)
    if scale is not None:
        scale = units(scale)
        scaling = (length_unit.to(scale) / scale).magnitude * scale
    else:
        scaling = length_unit
    return scaling
