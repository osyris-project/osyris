# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np


def _read_bound_key(infofile, ncpu):
    """
    Read hilbert bound key.
    """
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

    # return np.array(bound_key, dtype=float)
    return bound_key


def _btest(i, pos):
    return bool(i & (1 << pos))


def _hilbert3d(x, y, z, bit_length):

    # logical,dimension(0:3*bit_length-1)::i_bit_mask
    # logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
    # integer,dimension(0:7,0:1,0:11)::state_diagram

    # integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

    # if(bit_length>bit_size(bit_length))then
    #  write(*,*)'Maximum bit length=',bit_size(bit_length)
    #  write(*,*)'stop in hilbert3d'
    #  stop
    # endif

    i_bit_mask = np.zeros(3 * bit_length, dtype=bool)
    x_bit_mask = np.zeros(bit_length, dtype=bool)
    y_bit_mask = np.zeros_like(x_bit_mask)
    z_bit_mask = np.zeros_like(x_bit_mask)
    # logical,dimension(0:3*bit_length-1)::i_bit_mask
    # logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask

    # state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
    #    &   0, 1, 3, 2, 7, 6, 4, 5,&
    #    &   2, 6, 0, 7, 8, 8, 0, 7,&
    #    &   0, 7, 1, 6, 3, 4, 2, 5,&
    #    &   0, 9,10, 9, 1, 1,11,11,&
    #    &   0, 3, 7, 4, 1, 2, 6, 5,&
    #    &   6, 0, 6,11, 9, 0, 9, 8,&
    #    &   2, 3, 1, 0, 5, 4, 6, 7,&
    #    &  11,11, 0, 7, 5, 9, 0, 7,&
    #    &   4, 3, 5, 2, 7, 0, 6, 1,&
    #    &   4, 4, 8, 8, 0, 6,10, 6,&
    #    &   6, 5, 1, 2, 7, 4, 0, 3,&
    #    &   5, 7, 5, 3, 1, 1,11,11,&
    #    &   4, 7, 3, 0, 5, 6, 2, 1,&
    #    &   6, 1, 6,10, 9, 4, 9,10,&
    #    &   6, 7, 5, 4, 1, 0, 2, 3,&
    #    &  10, 3, 1, 1,10, 3, 5, 9,&
    #    &   2, 5, 3, 4, 1, 6, 0, 7,&
    #    &   4, 4, 8, 8, 2, 7, 2, 3,&
    #    &   2, 1, 5, 6, 3, 0, 4, 7,&
    #    &   7, 2,11, 2, 7, 5, 8, 5,&
    #    &   4, 5, 7, 6, 3, 2, 0, 1,&
    #    &  10, 3, 2, 6,10, 3, 4, 4,&
    #    &   6, 1, 7, 0, 5, 2, 4, 3 /), &
    #    & (/8 ,2, 12 /) )

    state_diagram = np.array([
        1, 2, 3, 2, 4, 5, 3, 5, 0, 1, 3, 2, 7, 6, 4, 5, 2, 6, 0, 7, 8, 8, 0, 7, 0, 7, 1,
        6, 3, 4, 2, 5, 0, 9, 10, 9, 1, 1, 11, 11, 0, 3, 7, 4, 1, 2, 6, 5, 6, 0, 6, 11,
        9, 0, 9, 8, 2, 3, 1, 0, 5, 4, 6, 7, 11, 11, 0, 7, 5, 9, 0, 7, 4, 3, 5, 2, 7, 0,
        6, 1, 4, 4, 8, 8, 0, 6, 10, 6, 6, 5, 1, 2, 7, 4, 0, 3, 5, 7, 5, 3, 1, 1, 11, 11,
        4, 7, 3, 0, 5, 6, 2, 1, 6, 1, 6, 10, 9, 4, 9, 10, 6, 7, 5, 4, 1, 0, 2, 3, 10, 3,
        1, 1, 10, 3, 5, 9, 2, 5, 3, 4, 1, 6, 0, 7, 4, 4, 8, 8, 2, 7, 2, 3, 2, 1, 5, 6,
        3, 0, 4, 7, 7, 2, 11, 2, 7, 5, 8, 5, 4, 5, 7, 6, 3, 2, 0, 1, 10, 3, 2, 6, 10, 3,
        4, 4, 6, 1, 7, 0, 5, 2, 4, 3
    ]).reshape(8, 2, 12)

    # convert to binary
    for i in range(bit_length):
        x_bit_mask[i] = _btest(x, i)
        y_bit_mask[i] = _btest(y, i)
        z_bit_mask[i] = _btest(z, i)

    # interleave bits
    for i in range(bit_length):
        i_bit_mask[3 * i + 2] = x_bit_mask[i]
        i_bit_mask[3 * i + 1] = y_bit_mask[i]
        i_bit_mask[3 * i] = z_bit_mask[i]

    # build Hilbert ordering using state diagram
    cstate = 0
    for i in range(bit_length - 1, 1, -1):
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

    # save Hilbert key as double precision real ??
    order = 0
    for i in range(3 * bit_length):
        b0 = 0
        if i_bit_mask[i]:
            b0 = 1
        order = order + b0 * (2**i)

    return order


def hilbert_cpu_list(bounding_box, lmax, levelmax, infofile, ncpu, ndim):

    # Read bound key
    bound_key = _read_bound_key(infofile=infofile, ncpu=ncpu)
    print(bound_key)

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
    #     imin = int(xmin * dble(maxdom))
    #     imax = imin + 1
    #     jmin = int(ymin * dble(maxdom))
    #     jmax = jmin + 1
    #     kmin = int(zmin * dble(maxdom))
    #     kmax = kmin + 1
    # endif

    dkey = (2**(levelmax + 1) / maxdom)**ndim
    ndom = 1
    if bit_length > 0:
        ndom = 8
    idom = [imin, imax] * 4
    jdom = [jmin, jmax] * 4
    kdom = [kmin, kmax] * 4

    # idom(1)=imin; idom(2)=imax
    # idom(3)=imin; idom(4)=imax
    # idom(5)=imin; idom(6)=imax
    # idom(7)=imin; idom(8)=imax
    # jdom(1)=jmin; jdom(2)=jmin
    # jdom(3)=jmax; jdom(4)=jmax
    # jdom(5)=jmin; jdom(6)=jmin
    # jdom(7)=jmax; jdom(8)=jmax
    # kdom(1)=kmin; kdom(2)=kmin
    # kdom(3)=kmin; kdom(4)=kmin
    # kdom(5)=kmax; kdom(6)=kmax
    # kdom(7)=kmax; kdom(8)=kmax

    # do i=1,ndom

    bounding_min = [0, 0, 0, 0, 0, 0, 0, 0]
    bounding_max = [0, 0, 0, 0, 0, 0, 0, 0]
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

    # ncpu_read=0
    cpu_list = []
    for i in range(ndom):
        for j in range(cpu_min[i], cpu_max[i] + 1):
            if j + 1 not in cpu_list:
                cpu_list.append(j + 1)

    print(cpu_list)
    return cpu_list


# else
#    ncpu_read=ncpu
#    do j=1,ncpu
#       cpu_list(j)=j
#    end do
# end  if