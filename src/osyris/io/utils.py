# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import glob
import struct
import numpy as np
from ..core import Array


def generate_fname(nout, path="", ftype="", cpuid=1, ext=""):

    if len(path) > 0:
        if path[-1] != "/":
            path = path + "/"

    if nout == -1:
        filelist = sorted(glob.glob(path + "output*"))
        number = filelist[-1].split("_")[-1]
    else:
        number = str(nout).zfill(5)

    infile = path + "output_" + number
    if len(ftype) > 0:
        infile += "/" + ftype + "_" + number
        if cpuid >= 0:
            infile += ".out" + str(cpuid).zfill(5)

    if len(ext) > 0:
        infile += ext

    return infile


def read_parameter_file(fname=None, delimiter="="):
    """
    Read info file and create dictionary
    """
    out = {}
    with open(fname) as f:
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


def make_vector_arrays(data):
    """
    Merge vector components in 2d arrays.
    """
    components = list("xyz"[:data.meta["ndim"]])
    if len(components) > 1:
        skip = []
        for key in list(data.keys()):
            # Make list of 3 components
            comp_list = None
            rawkey = None
            if key.endswith("_x") and key not in skip:
                rawkey = key[:-2]
                comp_list = [rawkey + "_" + c for c in components]
            if "_x_" in key and key not in skip:
                rawkey = key.replace("_x", "")
                comp_list = [key.replace("_x_", "_{}_".format(c)) for c in components]

            if comp_list is not None:
                if all([item in data for item in comp_list]):
                    data[rawkey] = Array(values=np.array(
                        [data[c].values for c in comp_list]).T,
                                         unit=data[key].unit)
                    for c in comp_list:
                        del data[c]
                        skip.append(c)


def find_max_amr_level(levelmax, select):
    """
    Test the selection function in `select` on the range of possible AMR levels
    to determine the max level to read.
    """
    possible_levels = np.arange(1, levelmax + 1, dtype=int)
    func_test = select["level"](possible_levels)
    inds = np.argwhere(func_test).ravel()
    return possible_levels[inds.max()]
