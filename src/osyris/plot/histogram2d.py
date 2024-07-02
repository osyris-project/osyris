# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)

from typing import Union

import numpy as np
from pint import Quantity

from ..core import Array, Layer, Plot
from ..core.tools import finmax, finmin, to_bin_centers
from .parser import get_norm, parse_layer
from .render import render


def _parse_limit(limit, x, logx, reduction):
    autox = False
    if limit is None:
        if reduction == "min":
            limit = finmin(x.values)
        elif reduction == "max":
            limit = finmax(x.values)
        else:
            raise RuntimeError(
                f"_parse_limit: unknown reduction operation {reduction}."
            )
        autox = True
    else:
        if isinstance(limit, Quantity):
            limit = limit.to(x.unit.units).magnitude
        if logx:
            limit = np.log10(limit)
    return limit, autox


def histogram2d(
    x: Array,
    y: Array,
    *layers,
    bins: Union[int, tuple] = 256,
    mode: str = None,
    logx: bool = False,
    logy: bool = False,
    loglog: bool = False,
    norm: str = None,
    filename: str = None,
    # resolution: Union[int, dict] = 256,
    operation: str = "sum",
    title: str = None,
    # xmin: float = None,
    # xmax: float = None,
    # ymin: float = None,
    # ymax: float = None,
    vmin: float = None,
    vmax: float = None,
    plot: bool = True,
    ax: object = None,
    **kwargs,
) -> Plot:
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
    xvals = x.values
    yvals = y.values

    if loglog:
        logx = logy = True

    xbins, ybins = (bins, bins) if isinstance(bins, int) else bins

    # Construct some bin edges
    edges = []
    for axis, binning in zip([xvals, yvals], [xbins, ybins]):
        if isinstance(binning, int):
            xmin = finmin(axis)
            xmax = finmax(axis)
            if logx:
                edges.append(
                    np.logspace(
                        np.log10(xmin),
                        np.nextafter(np.log10(xmax), np.inf),
                        binning + 1,
                    )
                )
            else:
                edges.append(np.linspace(xmin, np.nextafter(xmax, np.inf), binning + 1))
        else:
            edges.append(binning)

    xedges, yedges = edges

    # # if isinstance(xbins, int):
    # #     xmin = finmin(x.values)
    # #     xmax = finmax(x.values)
    # #     if logx:
    # #         xedges = np.logspace(np.log10(xmin), np.log10(xmax), xbins + 1)
    # #     else:
    # #         xedges = np.linspace(xmin, xmax, xbins + 1)
    # # else:
    # #     xedges = xbins
    # # if isinstance(xbins, int):
    # #     xmin = finmin(x.values)
    # #     xmax = finmax(x.values)
    # #     if logx:
    # #         xedges = np.logspace(np.log10(xmin), np.log10(xmax), xbins + 1)
    # #     else:
    # #         xedges = np.linspace(xmin, xmax, xbins + 1)
    # # else:
    # #     xedges = xbins

    # if logx:
    #     x = np.log10(x)
    # if logy:
    #     y = np.log10(y)

    # nx = resolution
    # ny = resolution

    # xmin, autoxmin = _parse_limit(xmin, x, logx, "min")
    # xmax, autoxmax = _parse_limit(xmax, x, logx, "max")
    # ymin, autoymin = _parse_limit(ymin, y, logy, "min")
    # ymax, autoymax = _parse_limit(ymax, y, logy, "max")

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

    # dx = xmax - xmin
    # dy = ymax - ymin
    # if autoxmin:
    #     xmin = xmin - 0.05 * dx
    # if autoxmax:
    #     xmax = xmax + 0.05 * dx
    # if autoymin:
    #     ymin = ymin - 0.05 * dy
    # if autoymax:
    #     ymax = ymax + 0.05 * dy

    # # Construct some bin edges and centers
    # if logx:
    #     xedges = np.logspace(xmin, xmax, nx + 1)
    # else:
    #     xedges = np.linspace(xmin, xmax, nx + 1)
    # if logy:
    #     yedges = np.logspace(ymin, ymax, ny + 1)
    # else:
    #     yedges = np.linspace(ymin, ymax, ny + 1)

    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    to_render = []
    to_process = []
    operations = []

    # If no layers are defined, make a layer for counting cells
    if len(layers) == 0:
        layers = [Layer(Array(values=np.ones_like(xvals), name="counts"))]

    # Digitize the x and y positions to get bin indices
    x_bin_indices = np.digitize(xvals, xedges) - 1
    y_bin_indices = np.digitize(yvals, yedges) - 1
    counts = np.zeros((len(ycenters), len(xcenters)))
    np.add.at(counts, (y_bin_indices, x_bin_indices), 1)
    mask = counts == 0

    for layer in layers:
        if isinstance(layer, Array):
            layer = Layer(layer)
        layer = parse_layer(
            layer,
            mode=mode,
            operation=operation,
            norm=norm,
            vmin=vmin,
            vmax=vmax,
            **kwargs,
        )
        layer.kwargs.update(
            norm=get_norm(norm=layer.norm, vmin=layer.vmin, vmax=layer.vmax)
        )
        # to_process.append(layer.data.norm.values)
        # operations.append(layer.operation)

        binned = np.zeros((len(ycenters), len(xcenters)))
        np.add.at(binned, (y_bin_indices, x_bin_indices), layer.data.norm.values)
        # np.max.at(binned, (y_bin_indices, x_bin_indices), layer.data.norm.values)
        # to_render[ind]["data"]

        if layer.operation == "mean":
            with np.errstate(invalid="ignore"):
                binned /= counts

        to_render.append(
            {
                "mode": layer.mode,
                "params": layer.kwargs,
                "unit": layer.data.unit,
                "name": layer.data.name,
                "data": np.ma.masked_where(mask, binned),
            }
        )

    # # Send to numba histogramming
    # binned, counts = hist2d(
    #     x=xvals,
    #     y=yvals,
    #     values=np.array(to_process),
    #     xmin=xmin,
    #     xmax=xmax,
    #     nx=nx,
    #     ymin=ymin,
    #     ymax=ymax,
    #     ny=ny,
    # )

    # mask = counts == 0
    # for ind in range(len(to_process)):
    #     if operations[ind] == "mean":
    #         with np.errstate(invalid="ignore"):
    #             binned[ind, ...] /= counts
    #     to_render[ind]["data"] = np.ma.masked_where(mask, binned[ind, ...])

    to_return = {
        "x": xcenters,
        "y": ycenters,
        "layers": to_render,
        "filename": filename,
    }

    if plot:
        figure = render(
            x=xcenters, y=ycenters, data=to_render, logx=logx, logy=logy, ax=ax
        )
        figure["ax"].set_xlabel(x.label)
        figure["ax"].set_ylabel(y.label)
        figure["ax"].set_title(title)
        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
