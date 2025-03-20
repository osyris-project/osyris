# SPDX-License-Identifier: BSD-3-Clause

from typing import Iterable, Literal, Union

import numpy as np

from ..core import Array, Layer, Plot
from ..core.tools import finmax, finmin, to_bin_centers
from .parser import get_norm, parse_layer
from .render import render


def _padded_limits(x: np.ndarray, log: bool, pad: float = 0.05):
    if log:
        xmin, xmax = np.log10([x[0], x[-1]])
        padx = (xmax - xmin) * pad
        xmin, xmax = 10.0 ** (xmin - padx), 10.0 ** (xmax + padx)
    else:
        xmin, xmax = x[0], x[-1]
        padx = (x[-1] - x[0]) * pad
        xmin, xmax = x[0] - padx, x[-1] + padx
    return xmin, xmax


def hist2d(
    x: Array,
    y: Array,
    *layers,
    bins: Union[int, Iterable, tuple] = 256,
    mode: str = None,
    logx: bool = False,
    logy: bool = False,
    loglog: bool = False,
    norm: str = None,
    filename: str = None,
    operation: Literal["sum", "mean"] = "sum",
    title: str = None,
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

    :param bins: The number of bins to use. Default is 256. Can also be an array of
        bin edges. If a single integer is passed, the bin edges are determined
        automatically. Can also be a tuple of two integers or arrays,
        representing the number of bins in the x and y directions.

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

    :param operation: The operation to apply inside the bins of the histogram.
        Possible values are ``'sum'`` and ``'mean'``. Default is ``'sum'``.

    :param title: The title of the figure. Default is ``None``.

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

    # Construct bin edges
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

    xcenters = to_bin_centers(xedges)
    ycenters = to_bin_centers(yedges)

    # If no layers are defined, make a layer for counting cells
    if len(layers) == 0:
        layers = [Layer(Array(values=np.ones_like(xvals), name="counts"))]

    # Digitize the x and y positions to get bin indices
    x_bin_indices = np.digitize(xvals, xedges) - 1
    y_bin_indices = np.digitize(yvals, yedges) - 1
    sel = x_bin_indices >= 0
    sel &= x_bin_indices < len(xedges) - 1
    sel &= y_bin_indices >= 0
    sel &= y_bin_indices < len(yedges) - 1
    x_bin_indices = x_bin_indices[sel]
    y_bin_indices = y_bin_indices[sel]

    counts = np.zeros((len(ycenters), len(xcenters)))
    np.add.at(counts, (y_bin_indices, x_bin_indices), 1)
    mask = counts == 0

    to_render = []
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

        binned = np.zeros((len(ycenters), len(xcenters)))
        np.add.at(binned, (y_bin_indices, x_bin_indices), layer.data.norm.values[sel])

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

        figure["ax"].set_xlim(*_padded_limits(xedges, logx))
        figure["ax"].set_ylim(*_padded_limits(yedges, logy))

        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
