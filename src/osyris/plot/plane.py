# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from pint.quantity import Quantity
from .slice import get_slice_direction
from .render import render
from .parser import parse_layer
from ..core import Plot, Array
from ..core.tools import to_bin_centers
from scipy.stats import binned_statistic_2d
from scipy import ndimage as nd


def plane(*layers,
          direction="z",
          dx=0.0,
          dy=0.0,
          fname=None,
          title=None,
          sinks=True,
          plot=True,
          mode=None,
          norm=None,
          vmin=None,
          vmax=None,
          operation="mean",
          origin=None,
          resolution=256,
          ax=None,
          **kwargs):
    """
    Plot a 2D slice through the data domain.
    """

    if isinstance(layers, Array):
        layers = [layers]

    to_process = []
    to_render = []
    operations = []
    for layer in layers:
        data, settings, params = parse_layer(layer,
                                             mode=mode,
                                             norm=norm,
                                             vmin=vmin,
                                             vmax=vmax,
                                             operation=operation,
                                             **kwargs)

        to_process.append(data)
        to_render.append({
            "mode": settings["mode"],
            "params": params,
            "unit": data.unit.units,
            "name": data.name
        })
        operations.append(settings["operation"])

    dataset = to_process[0].parent

    # Set window size
    if dy == 0.0:
        dy = dx
    if not isinstance(dx, Quantity):
        dx *= dataset["xyz"].unit
    if not isinstance(dy, Quantity):
        dy *= dataset["xyz"].unit

    dir_vecs, origin = get_slice_direction(direction=direction,
                                           dataset=dataset,
                                           dx=0.5 * (dx + dy),
                                           origin=origin)

    # Distance to the plane
    xyz = dataset["xyz"] - origin
    diagonal = dataset["dx"] * np.sqrt(dataset.meta["ndim"]) * 0.5
    dist1 = np.sum(xyz * dir_vecs[0],
                   axis=1)  # / np.linalg.norm(dir_vecs[0][1])

    # Distance from center
    dist2 = xyz - diagonal

    # Select only the cells in contact with the slice,
    # at a distance less than dx/2
    cube = np.ravel(
        np.where(
            np.logical_and(
                np.abs(dist1) <= 1.0001 * diagonal,
                np.abs(dist2.norm) <=
                max(dx.magnitude, dy.magnitude) * 0.5 * np.sqrt(2.0))))

    ncells = len(cube)

    # Project coordinates onto the plane by taking dot product with axes
    # vectors
    coords = xyz[cube]
    datax = np.inner(coords, dir_vecs[1])
    datay = np.inner(coords, dir_vecs[2])
    datadx = diagonal[cube] / np.sqrt(dataset.meta["ndim"])

    # Define slice extent and resolution
    xmin = -0.5 * dx.magnitude
    xmax = xmin + dx.magnitude
    ymin = -0.5 * dy.magnitude
    ymax = ymin + dy.magnitude
    nx = resolution
    ny = resolution
    # dpx = (xmax - xmin) / float(nx)
    # dpy = (ymax - ymin) / float(ny)
    # x = np.linspace(xmin + 0.5 * dpx, xmax - 0.5 * dpx, nx).magnitude
    # y = np.linspace(ymin + 0.5 * dpy, ymax - 0.5 * dpy, ny).magnitude

    # counts = np.zeros([ny, nx])
    scalar_layer = []
    to_binning = []
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream"]:
            if to_process[ind].ndim < 3:
                uv = to_process[ind].array[cube]
                # to_binning.append(uv[:, 0])
                # to_binning.append(uv[:, 1])
            else:
                uv = np.inner(to_process[ind].array.take(cube, axis=0),
                              dir_vecs[1:])
            w = None
            if "color" in to_render[ind]["params"]:
                if isinstance(to_render[ind]["params"]["color"], Array):
                    w = to_render[ind]["params"]["color"].norm
                elif isinstance(to_render[ind]["params"]["color"], np.ndarray):
                    w = to_render[ind]["params"]["color"]
            if w is None:
                w = np.linalg.norm(uv, axis=1)
            else:
                w = w.take(cube, axis=0)
            print("uv.shape", uv.shape)
            to_binning.append(uv[:, 0])
            to_binning.append(uv[:, 1])
            to_binning.append(w)
            # to_process[ind] = np.concatenate((uv, w.reshape(ncells, 1)),
            #                                  axis=1)

            scalar_layer.append(False)
            # to_render[ind]["data"] = np.zeros(
            #     [ny, nx, to_process[ind].shape[1]])
        else:
            # to_process[ind] = to_process[ind].array[cube]
            to_binning.append(to_process[ind].array[cube])
            # to_render[ind]["data"] = np.zeros([ny, nx])
            scalar_layer.append(True)

    # Buffer for counts
    to_binning.append(np.ones_like(to_binning[0]))

    # Construct some bin edges
    xedges = np.linspace(xmin, xmax, nx + 1)
    yedges = np.linspace(ymin, ymax, ny + 1)

    # In the contour plots, x and y are the centers of the cells, instead of
    # the edges.
    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    print(datax.shape)
    print(datay.shape)
    print(to_binning[0].shape)
    print(to_binning[0])
    # print(to_process[1].shape)
    print(xedges.shape)
    print(yedges.shape)
    print(len(to_binning))
    print(datax)
    print(xedges)
    print(type(xedges))
    print(type(yedges))

    # binned = []
    # for i in range(len(to_process)):
    #     binned[i], _, _ = np.histogram2d(datay,
    #                                      datax,
    #                                      bins=(yedges, xedges),
    #                                      weights=to_process[i])

    binned, _, _, _ = binned_statistic_2d(x=datay.values,
                                          y=datax.values,
                                          values=to_binning,
                                          statistic="mean",
                                          bins=[yedges, xedges])

    # print(binned[0])
    # mask = binned[0] == 0.0
    # print(mask)
    transform = tuple(
        nd.distance_transform_edt(np.isnan(binned[-1]),
                                  return_distances=False,
                                  return_indices=True))
    # inds = tuple(inds)
    counter = 0
    for ind in range(len(to_render)):
        # if operations[ind] != operation:
        #     if operation == "sum":
        #         with np.errstate(invalid="ignore"):
        #             binned[ind] /= binned[0]
        #     else:
        #         binned[ind] *= counts
        # binned[ind][mask] = np.nan
        # print(binned[ind])
        # inds = nd.distance_transform_edt(mask,
        #                                  return_distances=False,
        #                                  return_indices=True)

        # to_render[ind - 1]["data"] = np.ma.masked_where(mask, binned[ind])
        # print(inds)
        if scalar_layer[ind]:
            to_render[ind]["data"] = binned[counter][transform]
            counter += 1
        else:
            to_render[ind]["data"] = np.array([
                binned[counter + 1][transform], binned[counter][transform],
                binned[counter + 2][transform]
            ]).T
            counter += 3

        # print(to_render[ind - 1]["data"].shape)

    # ind = nd.distance_transform_edt(invalid,
    #                                 return_distances=False,
    #                                 return_indices=True)
    # return data[tuple(ind)]

    # datax -= xmin
    # datay -= ymin
    # istart = ((datax - datadx) / dpx).array.astype(np.int64)
    # iend = ((datax + datadx) / dpx).array.astype(np.int64) + 1
    # jstart = ((datay - datadx) / dpy).array.astype(np.int64)
    # jend = ((datay + datadx) / dpy).array.astype(np.int64) + 1

    # for i in range(len(istart)):
    #     i0 = istart[i]
    #     i1 = iend[i]
    #     j0 = jstart[i]
    #     j1 = jend[i]
    #     if i0 <= nx and j0 <= ny and i1 > 0 and j1 > 0:
    #         i0 = max(i0, 0)
    #         i1 = min(i1, nx)
    #         j0 = max(j0, 0)
    #         j1 = min(j1, ny)
    #         for ind in range(len(to_process)):
    #             to_render[ind]["data"][j0:j1, i0:i1] += to_process[ind][i]
    #         counts[j0:j1, i0:i1] += 1.0

    # # Normalize by counts
    # mask = counts == 0.0
    # for ind in range(len(to_process)):
    #     if to_render[ind]["data"].ndim > counts.ndim:
    #         to_render[ind]["data"] = np.ma.masked_where(
    #             np.broadcast_to(mask.reshape(ny, nx, 1),
    #                             to_render[ind]["data"].shape),
    #             to_render[ind]["data"]) / counts.reshape(ny, nx, 1)
    #     else:
    #         to_render[ind]["data"] = np.ma.masked_where(
    #             mask, to_render[ind]["data"]) / counts

    # import matplotlib.pyplot as plt
    # fig1, ax1 = plt.subplots()
    # ax1.imshow(counts, origin="lower", cmap="jet")
    # fig1.show()

    to_return = {"x": xcenters, "y": ycenters, "layers": to_render}
    if plot:
        # Render the map
        figure = render(x=xcenters, y=ycenters, data=to_render, ax=ax)
        figure["ax"].set_xlabel(dataset["xyz"].x.label)
        figure["ax"].set_ylabel(dataset["xyz"].y.label)
        if ax is None:
            figure["ax"].set_aspect("equal")
        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
