# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..core import Plot
from .. import units
from .render import render
from ..core.tools import to_bin_centers, finmin, finmax
from .parser import parse_layer
from scipy.stats import binned_statistic_2d


def histogram1d(x,
                weights=None,
                ax=None,
                logx=False,
                logy=False,
                loglog=False,
                norm=None,
                filename=None,
                bins=50,
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

    xvals = x.norm.values

    # Define plotting range
    autoxmin = False
    autoxmax = False
    autoymin = False
    autoymax = False

    if xmin is None:
        xmin = finmin(xvals)
        autoxmin = True
    if xmax is None:
        xmax = finmax(xvals)
        autoxmax = True

    if logx:
        [xmin, xmax] = np.log10([xmin, xmax])
    # if logy:
    #     [ymin, ymax] = np.log10([ymin, ymax])

    # # Protect against empty plots if xmin==xmax or ymin==ymax
    # if xmin == xmax:
    #     if xmin == 0.0:
    #         xmin = -0.1
    #         xmax = 0.1
    #     else:
    #         xmin = xmin - 0.05 * abs(xmin)
    #         xmax = xmax + 0.05 * abs(xmax)
    # if ymin == ymax:
    #     if ymin == 0.0:
    #         ymin = -0.1
    #         ymax = 0.1
    #     else:
    #         ymin = ymin - 0.05 * abs(ymin)
    #         ymax = ymax + 0.05 * abs(ymax)

    dx = xmax - xmin
    # dy = ymax - ymin
    if autoxmin:
        xmin = xmin - 0.05 * dx
    if autoxmax:
        xmax = xmax + 0.05 * dx
    # if autoymin:
    #     ymin = ymin - 0.05 * dy
    # if autoymax:
    #     ymax = ymax + 0.05 * dy

    # Construct some bin edges
    if logx:
        xedges = np.logspace(xmin, xmax, nx + 1)
    else:
        xedges = np.linspace(xmin, xmax, nx + 1)
    # if logy:
    #     yedges = np.logspace(ymin, ymax, ny + 1)
    # else:
    #     yedges = np.linspace(ymin, ymax, ny + 1)

    hist, _ = np.histogram1d(xvals, bins=xedges, weights=weights)
    if operation == "mean":
        counts, _ = np.histogram1d(np.ones_like(xvals), bins=xedges)
        hist /= counts

    # Here we assume that dictionary retains order of insertion: counts
    # are the first key
    mask = binned[0] == 0.0
    for ind in range(1, len(to_process)):
        if operations[ind] != operation:
            if operation == "sum":
                with np.errstate(invalid="ignore"):
                    binned[ind] /= binned[0]
            else:
                binned[ind] *= counts
        to_render[ind - 1]["data"] = np.ma.masked_where(mask, binned[ind])

    if len(to_render) == 0:
        _, _, params = parse_layer(layer=None,
                                   mode=mode,
                                   norm=norm,
                                   vmin=vmin,
                                   vmax=vmax,
                                   **kwargs)
        to_render.append({
            "data": np.ma.masked_where(mask, binned[0]),
            "mode": mode,
            "params": params,
            "unit": units.dimensionless,
            "name": "counts"
        })

    figure = render(x=xcenters, y=ycenters, data=to_render, logx=logx, logy=logy, ax=ax)

    figure["ax"].set_xlabel(x.label)
    figure["ax"].set_ylabel(y.label)

    return Plot(x=xcenters,
                y=ycenters,
                layers=to_render,
                fig=figure["fig"],
                ax=figure["ax"])
