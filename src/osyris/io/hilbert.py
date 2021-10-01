# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..core import Array


def _read_bound_key(infofile, ncpu):
    with open(infofile) as f:
        content = f.readlines()
    starting_line = None
    for num, line in enumerate(content):
        if len(set(["DOMAIN", "ind_min", "ind_max"]) - set(line.split())) == 0:
            starting_line = num + 1
            break
    bound_key = []
    if starting_line is not None:
        for n in range(ncpu):
            line = content[starting_line + n].split()
            bound_key.append(int(float(line[1])))
        bound_key.append(int(float(content[starting_line + ncpu - 1].split()[2])))
    return bound_key


def _btest(i, pos):
    return bool(i & (1 << pos))


def _hilbert3d(x, y, z, bit_length):

    i_bit_mask = np.zeros(3 * bit_length, dtype=bool)
    x_bit_mask = np.zeros(bit_length, dtype=bool)
    y_bit_mask = np.zeros_like(x_bit_mask)
    z_bit_mask = np.zeros_like(x_bit_mask)

    state_diagram = np.array([
        1, 2, 3, 2, 4, 5, 3, 5, 0, 1, 3, 2, 7, 6, 4, 5, 2, 6, 0, 7, 8, 8, 0, 7, 0, 7, 1,
        6, 3, 4, 2, 5, 0, 9, 10, 9, 1, 1, 11, 11, 0, 3, 7, 4, 1, 2, 6, 5, 6, 0, 6, 11,
        9, 0, 9, 8, 2, 3, 1, 0, 5, 4, 6, 7, 11, 11, 0, 7, 5, 9, 0, 7, 4, 3, 5, 2, 7, 0,
        6, 1, 4, 4, 8, 8, 0, 6, 10, 6, 6, 5, 1, 2, 7, 4, 0, 3, 5, 7, 5, 3, 1, 1, 11, 11,
        4, 7, 3, 0, 5, 6, 2, 1, 6, 1, 6, 10, 9, 4, 9, 10, 6, 7, 5, 4, 1, 0, 2, 3, 10, 3,
        1, 1, 10, 3, 5, 9, 2, 5, 3, 4, 1, 6, 0, 7, 4, 4, 8, 8, 2, 7, 2, 3, 2, 1, 5, 6,
        3, 0, 4, 7, 7, 2, 11, 2, 7, 5, 8, 5, 4, 5, 7, 6, 3, 2, 0, 1, 10, 3, 2, 6, 10, 3,
        4, 4, 6, 1, 7, 0, 5, 2, 4, 3
    ]).reshape((8, 2, 12), order='F')

    # Convert to binary
    for i in range(bit_length):
        x_bit_mask[i] = _btest(x, i)
        y_bit_mask[i] = _btest(y, i)
        z_bit_mask[i] = _btest(z, i)

    # Interleave bits
    for i in range(bit_length):
        i_bit_mask[3 * i + 2] = x_bit_mask[i]
        i_bit_mask[3 * i + 1] = y_bit_mask[i]
        i_bit_mask[3 * i] = z_bit_mask[i]

    # Build Hilbert ordering using state diagram
    cstate = 0
    for i in range(bit_length - 1, -1, -1):
        b2 = 0
        if i_bit_mask[3 * i + 2]:
            b2 = 1
        b1 = 0
        if i_bit_mask[3 * i + 1]:
            b1 = 1
        b0 = 0
        if i_bit_mask[3 * i]:
            b0 = 1
        sdigit = b2 * 4 + b1 * 2 + b0
        nstate = state_diagram[sdigit, 0, cstate]
        hdigit = state_diagram[sdigit, 1, cstate]
        i_bit_mask[3 * i + 2] = _btest(hdigit, 2)
        i_bit_mask[3 * i + 1] = _btest(hdigit, 1)
        i_bit_mask[3 * i] = _btest(hdigit, 0)
        cstate = nstate

    order = 0
    for i in range(3 * bit_length):
        b0 = 0
        if i_bit_mask[i]:
            b0 = 1
        order = order + b0 * (2**i)
    return order


def _get_cpu_list(bounding_box, lmax, levelmax, infofile, ncpu, ndim):

    bound_key = _read_bound_key(infofile=infofile, ncpu=ncpu)

    xmin = bounding_box["xmin"]
    xmax = bounding_box["xmax"]
    ymin = bounding_box["ymin"]
    ymax = bounding_box["ymax"]
    zmin = bounding_box["zmin"]
    zmax = bounding_box["zmax"]
    dmax = max(xmax - xmin, ymax - ymin, zmax - zmin)
    for ilevel in range(1, lmax + 1):
        dx = 0.5**ilevel
        if dx < dmax:
            break

    lmin = ilevel
    bit_length = lmin - 1
    maxdom = 2**bit_length
    imin = 0
    imax = 0
    jmin = 0
    jmax = 0
    kmin = 0
    kmax = 0
    if bit_length > 0:
        imin = int(xmin * maxdom)
        imax = imin + 1
        jmin = int(ymin * maxdom)
        jmax = jmin + 1
        kmin = int(zmin * maxdom)
        kmax = kmin + 1

    dkey = (2**(levelmax + 1) // maxdom)**ndim
    ndom = 1
    if bit_length > 0:
        ndom = 8
    idom = [imin, imax] * 4
    jdom = [jmin, jmin, jmax, jmax] * 2
    kdom = [kmin] * 4 + [kmax] * 4

    bounding_min = [0, 0, 0, 0, 0, 0, 0, 0]
    bounding_max = [0, 0, 0, 0, 0, 0, 0, 0]
    bounding = None
    for i in range(ndom):
        if bit_length > 0:
            bounding = _hilbert3d(idom[i], jdom[i], kdom[i], bit_length)
            order_min = bounding
        else:
            order_min = 0
        bounding_min[i] = order_min * dkey
        bounding_max[i] = (order_min + 1) * dkey

    cpu_min = [0, 0, 0, 0, 0, 0, 0, 0]
    cpu_max = [0, 0, 0, 0, 0, 0, 0, 0]
    for impi in range(ncpu):
        for i in range(ndom):
            if ((bound_key[impi] <= bounding_min[i])
                    and (bound_key[impi + 1] > bounding_min[i])):
                cpu_min[i] = impi

            if ((bound_key[impi] < bounding_max[i])
                    and (bound_key[impi + 1] >= bounding_max[i])):
                cpu_max[i] = impi

    cpu_list = []
    for i in range(ndom):
        for j in range(cpu_min[i], cpu_max[i] + 1):
            if j + 1 not in cpu_list:
                cpu_list.append(j + 1)

    return cpu_list


def hilbert_cpu_list(meta, scaling, select, infofile):
    if meta["ordering type"] != "hilbert":
        return
    bounding_box = {"xmin": 0, "xmax": 1, "ymin": 0, "ymax": 1, "zmin": 0, "zmax": 1}
    # Make an array of cell centers according to lmax
    box_size = (meta["boxlen"] * scaling).magnitude
    ncells = 2**min(meta["levelmax"], 18)  # limit to 262000 cells
    half_dxmin = 0.5 * box_size / ncells
    xyz_centers = Array(values=np.linspace(half_dxmin, box_size - half_dxmin, ncells),
                        unit=1.0 * scaling.units)
    new_bbox = False
    for c in "xyz":
        if c in select:
            new_bbox = True
            func_test = select[c](xyz_centers)
            inds = np.argwhere(func_test).ravel()
            start = xyz_centers[inds.min()] - (half_dxmin * scaling.units)
            end = xyz_centers[inds.max()] + (half_dxmin * scaling.units)
            bounding_box["{}min".format(c)] = start._array / box_size
            bounding_box["{}max".format(c)] = end._array / box_size
            select["xyz_{}".format(c)] = select.pop(c)

    if new_bbox:
        return _get_cpu_list(bounding_box=bounding_box,
                             lmax=meta["lmax"],
                             levelmax=meta["levelmax"],
                             infofile=infofile,
                             ncpu=meta["ncpu"],
                             ndim=meta["ndim"])
