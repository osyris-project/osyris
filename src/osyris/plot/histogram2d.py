# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..core import Array, Plot
from .render import render
from ..core.tools import to_bin_centers, finmin, finmax
from .parser import parse_layer
from .utils import hist2d


def histogram2d(x,
                y,
                *layers,
                mode=None,
                ax=None,
                logx=False,
                logy=False,
                loglog=False,
                norm=None,
                filename=None,
                resolution=256,
                operation="sum",
                title=None,
                xmin=None,
                xmax=None,
                ymin=None,
                ymax=None,
                vmin=None,
                vmax=None,
                **kwargs):
    """
    Plot a 2D histogram with two variables as input.
    """
    if loglog:
        logx = logy = True

    nx = resolution
    ny = resolution

    xvals = x.norm.values
    yvals = y.norm.values
    if logx:
        xvals = np.log10(xvals)
    if logy:
        yvals = np.log10(yvals)

    # Define plotting range
    autoxmin = False
    autoxmax = False
    autoymin = False
    autoymax = False

    if xmin is None:
        xmin = finmin(xvals)
        autoxmin = True
    else:
        xmin = np.log10(xmin)
    if xmax is None:
        xmax = finmax(xvals)
        autoxmax = True
    else:
        xmax = np.log10(xmax)
    if ymin is None:
        ymin = finmin(yvals)
        autoymin = True
    else:
        ymin = np.log10(ymin)
    if ymax is None:
        ymax = finmax(yvals)
        autoymax = True
    else:
        ymax = np.log10(ymax)

    # Protect against empty plots if xmin==xmax or ymin==ymax
    if xmin == xmax:
        if xmin == 0.0:
            xmin = -0.1
            xmax = 0.1
        else:
            xmin = xmin - 0.05 * abs(xmin)
            xmax = xmax + 0.05 * abs(xmax)
    if ymin == ymax:
        if ymin == 0.0:
            ymin = -0.1
            ymax = 0.1
        else:
            ymin = ymin - 0.05 * abs(ymin)
            ymax = ymax + 0.05 * abs(ymax)

    dx = xmax - xmin
    dy = ymax - ymin
    if autoxmin:
        xmin = xmin - 0.05 * dx
    if autoxmax:
        xmax = xmax + 0.05 * dx
    if autoymin:
        ymin = ymin - 0.05 * dy
    if autoymax:
        ymax = ymax + 0.05 * dy

    # Construct some bin edges and centers
    if logx:
        xedges = np.logspace(xmin, xmax, nx + 1)
    else:
        xedges = np.linspace(xmin, xmax, nx + 1)
    if logy:
        yedges = np.logspace(ymin, ymax, ny + 1)
    else:
        yedges = np.linspace(ymin, ymax, ny + 1)
    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    to_render = []
    to_process = []
    operations = []

    # If no layers are defined, make a layer for counting cells
    if len(layers) == 0:
        layers = [Array(values=np.ones_like(xvals), name="counts")]

    for layer in layers:
        data, settings, params = parse_layer(layer=layer,
                                             mode=mode,
                                             norm=norm,
                                             vmin=vmin,
                                             vmax=vmax,
                                             operation=operation,
                                             **kwargs)
        to_process.append(data.norm.values)
        to_render.append({
            "mode": settings["mode"],
            "params": params,
            "unit": data.unit.units,
            "name": data.name
        })
        operations.append(settings["operation"])

    # Send to numba histogramming
    binned, counts = hist2d(x=xvals,
                            y=yvals,
                            values=np.array(to_process),
                            xmin=xmin,
                            xmax=xmax,
                            nx=nx,
                            ymin=ymin,
                            ymax=ymax,
                            ny=ny)

    mask = counts == 0
    for ind in range(len(to_process)):
        if operations[ind] == "mean":
            with np.errstate(invalid="ignore"):
                binned[ind, ...] /= counts
        to_render[ind]["data"] = np.ma.masked_where(mask, binned[ind, ...])

    figure = render(x=xcenters, y=ycenters, data=to_render, logx=logx, logy=logy, ax=ax)

    figure["ax"].set_xlabel(x.label)
    figure["ax"].set_ylabel(y.label)

    return Plot(x=xcenters,
                y=ycenters,
                layers=to_render,
                fig=figure["fig"],
                ax=figure["ax"])
