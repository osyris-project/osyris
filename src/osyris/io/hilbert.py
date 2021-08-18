# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np


def hilbert_cpu_list(bounding_box, lmax, levelmax, infofile, ncpu):

    # Read bound key
    bound_key = read_bound_key(infofile=infofile, ncpu=ncpu)

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
    idom = [[imin, imax]] * 4
    jdom = [[jmin, jmax]] * 4
    kdom = [[kmin, kmax]] * 4

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


#    do i=1,ndom
#       if(bit_length>0)then
#          call hilbert3d(idom(i),jdom(i),kdom(i),bounding(1),bit_length,1)
#          order_min=bounding(1)
#       else
#          order_min=0.0d0
#       endif
#       bounding_min(i)=(order_min)*dkey
#       bounding_max(i)=(order_min+1.0D0)*dkey
#    end do

#    cpu_min=0; cpu_max=0
#    do impi=1,ncpu
#       do i=1,ndom
#          if (   bound_key(impi-1).le.bounding_min(i).and.&
#               & bound_key(impi  ).gt.bounding_min(i))then
#             cpu_min(i)=impi
#          endif
#          if (   bound_key(impi-1).lt.bounding_max(i).and.&
#               & bound_key(impi  ).ge.bounding_max(i))then
#             cpu_max(i)=impi
#          endif
#       end do
#    end do

#    ncpu_read=0
#    do i=1,ndom
#       do j=cpu_min(i),cpu_max(i)
#          if(.not. cpu_read(j))then
#             ncpu_read=ncpu_read+1
#             cpu_list(ncpu_read)=j
#             cpu_read(j)=.true.
#          endif
#       enddo
#    enddo
# else
#    ncpu_read=ncpu
#    do j=1,ncpu
#       cpu_list(j)=j
#    end do
# end  if


def read_bound_key(infofile, ncpu):
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
            bound_key.append(line[1])
        bound_key.append(content[starting_line + ncpu - 1].split()[2])

    return np.array(bound_key, dtype=float)
