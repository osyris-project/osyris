# SPDX-License-Identifier: BSD-3-Clause

from typing import Iterable, Union

import numpy as np

from ..core import Array, Layer, Plot
from ..core.tools import finmax, finmin, to_bin_centers
from .parser import parse_layer
from .render import render


def hist1d(
    *layers: Union[Iterable, Array],
    bins: Union[int, Iterable] = 50,
    weights: Array = None,
    logx: bool = False,
    logy: bool = False,
    loglog: bool = False,
    filename: str = None,
    title: str = None,
    ymin: float = None,
    ymax: float = None,
    ax: object = None,
    **kwargs,
) -> Plot:
    """
    Plot a 1D histogram with arbitrary number of variables as input.
    When a vector quantity is supplied, the function will histogram the norm of
    the vectors.


    This function has an API very close to that of matplotlib's ``hist`` function.
    For the documentation of any parameters that are not listed below, see
    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html.

    :param layers: Dicts or Arrays representing the quantities to be mapped onto the
        colormap of the generated image.

    :param bins: The number of bins to use. Default is 50. Can also be an array of
        bin edges. If a single integer is passed, the bin edges are determined
        automatically.

    :param weights: An array of the same length as the data, representing the weights
        of each data point. Default is ``None``, meaning all weights are 1.

    :param logx: If ``True``, use logarithmic scaling on the horizontal axis.
        Default is ``False``.

    :param logy: If ``True``, use logarithmic scaling on the vertical axis.
        Default is ``False``.

    :param loglog: If ``True``, use logarithmic scaling on the horizontal and
        vertical axes. Default is ``False``.

    :param filename: If specified, the returned figure is also saved to file.
        Default is ``None``.

    :param title: The title of the figure. Default is ``None``.

    :param ymin: Minimum value for the vertical axis. Default is ``None``.

    :param ymax: Maximum value for the vertical axis. Default is ``None``.

    :param ax: A matplotlib axes inside which the figure will be plotted.
        Default is ``None``, in which case some new axes a created.
    """
    if loglog:
        logx = logy = True

    figure = render(logx=logx, logy=logy, ax=ax)

    for layer in layers:
        if isinstance(layer, Array):
            layer = Layer(layer)
        layer = parse_layer(layer, bins=bins, weights=weights, **kwargs)

        xvals = layer.data.norm.values
        if layer.weights is not None:
            layer.weights = layer.weights.norm.values

        # Construct some bin edges
        if isinstance(layer.bins, int):
            xmin = finmin(xvals)
            xmax = finmax(xvals)
            if logx:
                xedges = np.logspace(
                    np.log10(xmin), np.nextafter(np.log10(xmax), np.inf), layer.bins + 1
                )
            else:
                xedges = np.linspace(xmin, np.nextafter(xmax, np.inf), layer.bins + 1)
        else:
            xedges = layer.bins

        ydata, _, _ = figure["ax"].hist(
            xvals, bins=xedges, weights=layer.weights, **layer.kwargs
        )

        figure["ax"].set_xlabel(layer.data.label)

    figure["ax"].set_ylim(ymin, ymax)
    figure["ax"].set_title(title)
    return Plot(
        x=to_bin_centers(xedges),
        y=ydata,
        fig=figure["fig"],
        ax=figure["ax"],
        filename=filename,
    )
