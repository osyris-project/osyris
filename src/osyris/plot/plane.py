# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
import numpy.ma as ma
from pint.quantity import Quantity
from .slice import get_slice_direction
from .render import render
from .parser import parse_layer
from ..core import Plot, Array
from ..core.tools import to_bin_centers, apply_mask
from scipy.stats import binned_statistic_2d
from scipy.ndimage import distance_transform_edt


def plane(*layers,
          direction="z",
          dx=None,
          dy=None,
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
    if dy is None:
        dy = dx
    if dx is not None and not isinstance(dx, Quantity):
        dx *= dataset["xyz"].unit
    if dy is not None and not isinstance(dy, Quantity):
        dy *= dataset["xyz"].unit

    dir_vecs, origin = get_slice_direction(direction=direction,
                                           dataset=dataset,
                                           dx=dx,
                                           dy=dy,
                                           origin=origin)

    # Distance to the plane
    xyz = dataset["xyz"] - origin
    diagonal = dataset["dx"] * np.sqrt(dataset.meta["ndim"]) * 0.5
    dist1 = np.sum(xyz * dir_vecs[0], axis=1)  # / np.linalg.norm(dir_vecs[0][1])

    # Select cells in contact with plane
    cells_in_plane = np.abs(dist1) <= diagonal
    # Project coordinates onto the plane by taking dot product with axes
    # vectors
    select = np.ravel(np.where(cells_in_plane))
    coords = xyz[select]
    datax = np.inner(coords, dir_vecs[1])
    datay = np.inner(coords, dir_vecs[2])
    datadx = 0.5 * dataset["dx"][select]

    # Get limits
    limits = {
        'xmin': np.amin(datax - datadx).values,
        'xmax': np.amax(datax + datadx).values,
        'ymin': np.amin(datay - datadx).values,
        'ymax': np.amax(datay + datadx).values
    }

    # Define slice extent
    if dx is None:
        xmin = limits['xmin']
        xmax = limits['xmax']
        ymin = limits['ymin']
        ymax = limits['ymax']
    else:
        xmin = -0.5 * dx.magnitude
        xmax = xmin + dx.magnitude
        ymin = -0.5 * dy.magnitude
        ymax = ymin + dy.magnitude
        # Limit selection further by using distance from center
        dist2 = xyz - diagonal
        select = np.ravel(
            np.where(
                np.logical_and(
                    cells_in_plane,
                    np.abs(dist2.norm) <=
                    max(dx.magnitude, dy.magnitude) * 0.5 * np.sqrt(2.0))))
        coords = xyz[select]
        datax = np.inner(coords, dir_vecs[1])
        datay = np.inner(coords, dir_vecs[2])

    scalar_layer = []
    to_binning = []
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream"]:
            if to_process[ind].ndim < 3:
                uv = to_process[ind].array[select]
            else:
                uv = np.inner(to_process[ind].array.take(select, axis=0), dir_vecs[1:])
            w = None
            if "color" in to_render[ind]["params"]:
                if isinstance(to_render[ind]["params"]["color"], Array):
                    w = to_render[ind]["params"]["color"].norm
                elif isinstance(to_render[ind]["params"]["color"], np.ndarray):
                    w = to_render[ind]["params"]["color"]
            if w is None:
                w = np.linalg.norm(uv, axis=1)
            else:
                w = w.take(select, axis=0)
            to_binning.append(apply_mask(uv[:, 0]))
            to_binning.append(apply_mask(uv[:, 1]))
            to_binning.append(apply_mask(w))
            scalar_layer.append(False)
        else:
            to_binning.append(apply_mask(to_process[ind].array[select]))
            scalar_layer.append(True)

    # Buffer for counts
    to_binning.append(np.ones_like(to_binning[0]))

    # Construct some bin edges
    xedges = np.linspace(xmin, xmax, resolution + 1)
    yedges = np.linspace(ymin, ymax, resolution + 1)

    # In the contour plots, x and y are the centers of the cells, instead of
    # the edges.
    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    binned, _, _, _ = binned_statistic_2d(x=apply_mask(datay.array),
                                          y=apply_mask(datax.array),
                                          values=to_binning,
                                          statistic="mean",
                                          bins=[yedges, xedges])

    # Use Scipy's distance transform to fill blanks with nearest neighbours
    condition = np.isnan(binned[-1])
    transform = tuple(
        distance_transform_edt(condition, return_distances=False, return_indices=True))

    # Define a mask to mask aread outside of the domain range, which have
    # been filled by the previous transform step
    xx = np.broadcast_to(xcenters, [resolution, resolution])
    yy = np.broadcast_to(ycenters, [resolution, resolution]).T
    mask = np.logical_or.reduce((xx < limits['xmin'], xx > limits['xmax'],
                                 yy < limits['ymin'], yy > limits['ymax']))
    mask_vec = np.broadcast_to(mask.reshape(*mask.shape, 1), mask.shape + (3, ))

    counter = 0
    for ind in range(len(to_render)):
        if scalar_layer[ind]:
            to_render[ind]["data"] = ma.masked_where(mask,
                                                     binned[counter][transform],
                                                     copy=False)
            counter += 1
        else:
            to_render[ind]["data"] = ma.masked_where(
                mask_vec,
                np.array([
                    binned[counter][transform].T, binned[counter + 1][transform].T,
                    binned[counter + 2][transform].T
                ]).T,
                copy=False)
            counter += 3

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
