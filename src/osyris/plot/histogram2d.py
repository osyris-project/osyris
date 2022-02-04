# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from pint.quantity import Quantity
from typing import Union
from ..core import Array, Plot
from .render import render
from ..core.tools import to_bin_centers, finmin, finmax
from .parser import parse_layer
from .utils import hist2d


def _parse_limit(limit, x, logx, reduction):
    autox = False
    if limit is None:
        if reduction == "min":
            limit = finmin(x.values)
        elif reduction == "max":
            limit = finmax(x.values)
        else:
            raise RuntimeError(
                f"_parse_limit: unknown reduction operation {reduction}.")
        autox = True
    else:
        if isinstance(limit, Quantity):
            limit = limit.to(x.unit.units).magnitude
        if logx:
            limit = np.log10(limit)
    return limit, autox


def histogram2d(x: Array,
                y: Array,
                *layers,
                mode: str = None,
                logx: bool = False,
                logy: bool = False,
                loglog: bool = False,
                norm: str = None,
                filename: str = None,
                resolution: Union[int, dict] = 256,
                operation: str = "sum",
                title: str = None,
                xmin: float = None,
                xmax: float = None,
                ymin: float = None,
                ymax: float = None,
                vmin: float = None,
                vmax: float = None,
                plot: bool = True,
                ax: object = None,
                **kwargs) -> Plot:
    """
    Plot a 2D histogram with two variables as input.
    When a vector quantity is supplied, the function will histogram the norm of
    the vectors.


    :param x: Horizontal Array to be histogrammed.

    :param y: Vertical Array to be histogrammed.

    :param layers: Dicts or Arrays representing the quantities to be mapped onto the
        colormap of the generated image.

    :param mode: The rendering mode for the histogram. Possible choices are
        ``'image'``, ``'contourf'``, ``'contour'``, and ``'scatter'``. Default is
        ``None``, which selects the ``render_mode`` set in the user configuration
        file (``'image'`` by default).

    :param logx: If ``True``, use logarithmic scaling on the horizontal axis.
        Default is ``False``.

    :param logy: If ``True``, use logarithmic scaling on the vertical axis.
        Default is ``False``.

    :param loglog: If ``True``, use logarithmic scaling on the horizontal and
        vertical axes. Default is ``False``.

    :param norm: The colormap normalization. Possible values are ``'linear'`` and
        ``'log'``. Default is ``None`` (= ``'linear'``).

    :param filename: If specified, the returned figure is also saved to file.
        Default is ``None``.

    :param resolution: Resolution of the generated map. This can either be an
        integer or a dict. In the case of an integer, it represents the number of
        pixels used for the horizontal and vertical dimensions. For a dictionary,
        the following syntax should be used: ``resolution={'x': 128, 'y': 192}``.
        Default is ``256``.

    :param operation: The operation to apply inside the bins of the histogram.
        Possible values are ``'sum'`` and ``'mean'``. Default is ``'sum'``.

    :param title: The title of the figure. Default is ``None``.

    :param xmin: Minimum value for the horizontal axis. Default is ``None``.

    :param xmax: Maximum value for the horizontal axis. Default is ``None``.

    :param ymin: Minimum value for the vertical axis. Default is ``None``.

    :param ymax: Maximum value for the vertical axis. Default is ``None``.

    :param vmin: Minimum value for colorbar range. Default is ``None``.

    :param vmax: Maximum value for colorbar range. Default is ``None``.

    :param plot: Make a plot if ``True``. If not, just return the ``Plot`` object
        containing the data that would be used to generate the plot.
        Default is ``True``.

    :param ax: A matplotlib axes inside which the figure will be plotted.
        Default is ``None``, in which case some new axes a created.
    """

    x = x.norm
    y = y.norm

    if loglog:
        logx = logy = True
    if logx:
        x = np.log10(x)
    if logy:
        y = np.log10(y)

    nx = resolution
    ny = resolution

    xmin, autoxmin = _parse_limit(xmin, x, logx, "min")
    xmax, autoxmax = _parse_limit(xmax, x, logx, "max")
    ymin, autoymin = _parse_limit(ymin, y, logy, "min")
    ymax, autoymax = _parse_limit(ymax, y, logy, "max")

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

    xvals = x.values
    yvals = y.values

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

    to_return = {
        "x": xcenters,
        "y": ycenters,
        "layers": to_render,
        "filename": filename
    }

    if plot:
        figure = render(x=xcenters,
                        y=ycenters,
                        data=to_render,
                        logx=logx,
                        logy=logy,
                        ax=ax)
        figure["ax"].set_xlabel(x.label)
        figure["ax"].set_ylabel(y.label)
        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
