# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from pint.quantity import Quantity
from .slice import get_slice_direction
from .render import render
from .parser import parse_layer
from ..core import Plot, Array


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
          origin=[0, 0, 0],
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
    diagonal = dataset["dx"] * np.sqrt(3.0) * 0.5
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
    datadx = diagonal[cube]

    # Define slice extent and resolution
    xmin = -0.5 * dx
    xmax = xmin + dx
    ymin = -0.5 * dy
    ymax = ymin + dy
    nx = resolution
    ny = resolution
    dpx = (xmax - xmin) / float(nx)
    dpy = (ymax - ymin) / float(ny)
    x = np.linspace(xmin + 0.5 * dpx, xmax - 0.5 * dpx, nx).magnitude
    y = np.linspace(ymin + 0.5 * dpy, ymax - 0.5 * dpy, ny).magnitude

    counts = np.zeros([ny, nx])
    for ind in range(len(to_process)):
        if to_render[ind]["mode"] in ["vec", "stream"]:
            if to_process[ind].ndim < 3:
                to_process[ind] = to_process[ind].array[cube]
            else:
                uv = np.inner(to_process[ind].array.take(cube, axis=0),
                              dir_vecs[1:])
                w = None
                if "color" in to_render[ind]["params"]:
                    if isinstance(to_render[ind]["params"]["color"], Array):
                        w = to_render[ind]["params"]["color"].norm
                    elif isinstance(to_render[ind]["params"]["color"],
                                    np.ndarray):
                        w = to_render[ind]["params"]["color"]
                if w is None:
                    w = np.linalg.norm(uv, axis=1)
                else:
                    w = w.take(cube, axis=0)
                to_process[ind] = np.concatenate((uv, w.reshape(ncells, 1)),
                                                 axis=1)
            to_render[ind]["data"] = np.zeros(
                [ny, nx, to_process[ind].shape[1]])
        else:
            to_process[ind] = to_process[ind].array[cube]
            to_render[ind]["data"] = np.zeros([ny, nx])

    datax -= xmin
    datay -= ymin
    istart = ((datax - datadx) / dpx).array.astype(np.int64)
    iend = ((datax + datadx) / dpx).array.astype(np.int64) + 1
    jstart = ((datay - datadx) / dpy).array.astype(np.int64)
    jend = ((datay + datadx) / dpy).array.astype(np.int64) + 1

    for i in range(len(istart)):
        i0 = istart[i]
        i1 = iend[i]
        j0 = jstart[i]
        j1 = jend[i]
        if i0 <= nx and j0 <= ny and i1 > 0 and j1 > 0:
            i0 = max(i0, 0)
            i1 = min(i1, nx)
            j0 = max(j0, 0)
            j1 = min(j1, ny)
            for ind in range(len(to_process)):
                to_render[ind]["data"][j0:j1, i0:i1] += to_process[ind][i]
            counts[j0:j1, i0:i1] += 1.0

    # Normalize by counts
    mask = counts == 0.0
    for ind in range(len(to_process)):
        if to_render[ind]["data"].ndim > counts.ndim:
            to_render[ind]["data"] = np.ma.masked_where(
                np.broadcast_to(mask.reshape(ny, nx, 1),
                                to_render[ind]["data"].shape),
                to_render[ind]["data"]) / counts.reshape(ny, nx, 1)
        else:
            to_render[ind]["data"] = np.ma.masked_where(
                mask, to_render[ind]["data"]) / counts

    to_return = {"x": x, "y": y, "layers": to_render}
    if plot:
        # Render the map
        figure = render(x=x, y=y, data=to_render, ax=ax)
        figure["ax"].set_xlabel(dataset["xyz"].x.label)
        figure["ax"].set_ylabel(dataset["xyz"].y.label)
        if ax is None:
            figure["ax"].set_aspect("equal")
        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
