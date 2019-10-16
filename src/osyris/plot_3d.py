# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
try:
    import ipyvolume as ipv
    no_3d = False
except ImportError:
    print("Warning: 3d plots are disabled because ipyvolume was not found.")
    no_3d = True


def plot_volume(field=False, dx=0.0, dy=0.0, dz=0.0, fname=None, title=None,
                sinks=True, resolution=128, **kwargs):
    """
    Plot volume redering of 3D data cube.
    """

    if no_3d:
        print("plot_volume is disabled because ipyvolume is not installed.")
        return

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


def plot_quiver(field=False, iskip=1, dx=0.0, dy=0.0, dz=0.0, fname=None,
                title=None, sinks=True, size=1, **kwargs):

    if no_3d:
        print("plot_quiver is disabled because ipyvolume is not installed.")
        return

    # Find parent container of object to plot
    holder = field.parent

    cube = np.where(np.logical_and(holder.get("x") >= -0.5*dx, \
                    np.logical_and(holder.get("x") <=  0.5*dx, \
                    np.logical_and(holder.get("y") >= -0.5*dy, \
                    np.logical_and(holder.get("y") <=  0.5*dy, \
                    np.logical_and(holder.get("z") >= -0.5*dz, \
                                   holder.get("z") <=  0.5*dz))))))

    ipv.figure()
    quiver = ipv.quiver(holder.get("x")[cube][::iskip],
                        holder.get("y")[cube][::iskip],
                        holder.get("z")[cube][::iskip],
                        field.x.get()[cube][::iskip],
                        field.y.get()[cube][::iskip],
                        field.z.get()[cube][::iskip], size=size)
    ipv.show()
    return
