# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np


def plot_volume(field=False, dx=0.0, dy=0.0, dz=0.0, fname=None, title=None,
                sinks=True, resolution=128, **kwargs):
    """
    Plot a 2D slice through the data cube. The arguments are:
    - scalar     : the scalar field to be plotted, e.g. mydata.density
    - image      : the scalar field to be plotted with an image
    - contour    : the scalar field to be plotted with contours
    - vec        : the vector field to be overplotted, e.g. mydata.velocity
    - stream     : the field for streamlines to be overplotted, e.g. mydata.B
    - direction  : the direction normal to the plane of the slice
    - dx         : the x extent of the slice, in units of scale (see data loader)
    - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
    - dz         : the thickness of the slice
    - axes       : if specified, the data is plotted on the specified axes (see demo).
    - resolution : number of pixels in the slice.
    - fname      : if specified, the figure is saved to file.
    """

    import ipyvolume as ipv

    # Find parent container of object to plot
    holder = field.parent

    cube = np.where(np.logical_and(holder.get("x") >= -0.5*dx, \
                    np.logical_and(holder.get("x") <=  0.5*dx, \
                    np.logical_and(holder.get("y") >= -0.5*dy, \
                    np.logical_and(holder.get("y") <=  0.5*dy, \
                    np.logical_and(holder.get("z") >= -0.5*dz, \
                                   holder.get("z") <=  0.5*dz))))))

    V0, edges = np.histogramdd(np.array([holder.get("x")[cube],
                                         holder.get("y")[cube],
                                         holder.get("z")[cube]]).T,
                               bins=(resolution, resolution, resolution))
    V1, edges = np.histogramdd(np.array([holder.get("x")[cube],
                                         holder.get("y")[cube],
                                         holder.get("z")[cube]]).T,
                               bins=(resolution, resolution, resolution),
                               weights=field.get()[cube])
    V = V1 / V0
    ipv.quickvolshow(V)
    ipv.show()

    return
