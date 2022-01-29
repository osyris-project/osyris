# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from pint.quantity import Quantity
from ..core import Plot, Array
from .. import units
from .render import render
from .parser import parse_layer


def line(x: Array,
         *y,
         logx: bool = False,
         logy: bool = False,
         loglog: bool = False,
         filename: str = None,
         title: str = None,
         xmin: float = None,
         xmax: float = None,
         ymin: float = None,
         ymax: float = None,
         ax: object = None,
         **kwargs) -> Plot:
    """
    Make a 1D plot with two variables as input.

    This function has an API very close to that of matplotlib's ``scatter`` function.
    For the documentation of any parameters that are not listed below, see
    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html.

    :param x: Array to use for scatter point positions along the horizontal dimension.

    :param y: Array to use for scatter point positions along the vertical dimension.

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

    :param title: The title of the figure. Default is ``None``.

    :param xmin: Minimum value for the horizontal axis. Default is ``None``.

    :param xmax: Maximum value for the horizontal axis. Default is ``None``.

    :param ymin: Minimum value for the vertical axis. Default is ``None``.

    :param ymax: Maximum value for the vertical axis. Default is ``None``.

    :param vmin: Minimum value for colorbar range. Default is ``None``.

    :param vmax: Maximum value for colorbar range. Default is ``None``.

    :param ax: A matplotlib axes inside which the figure will be plotted.
        Default is ``None``, in which case some new axes a created.
    """
    if loglog:
        logx = logy = True

    to_render = []
    yaxis_unit = None

    if isinstance(x, dict):
        xvals = x["x"].norm.values
        to_render.append({
            "data": x["y"].norm.values,
            "mode": "line",
            "params": {
                **kwargs
            },
            "unit": x["y"].unit.units,
            "name": x["y"].name
        })
        yaxis_unit = x["y"].unit
    else:
        xvals = x.norm.values

    sorting = np.argsort(xvals)

    if isinstance(layers, Array):
        layers = [layers]

    # to_process = []
    # to_render = []
    # to_scatter = []
    # yaxis_unit = None
    for layer in layers:
        data, settings, params = parse_layer(layer=layer, **kwargs)
        # to_process.append(data)
        del params["norm"]
        to_render.append({
            "data": data.norm.values[sorting],
            "mode": "line",
            "params": params,
            "unit": data.unit.units,
            "name": data.name
        })
        if yaxis_unit is None:
            yaxis_unit = data.unit
        else:
            if data.unit != yaxis_unit:
                raise RuntimeError(
                    "Different layers in 1D plots must all have the same unit.")

    figure = render(x=xvals[sorting], data=to_render, logx=logx, logy=logy, ax=ax)

    figure["ax"].set_xlabel(x.label)
    figure["ax"].set_ylabel(layers[0].label)
    figure["ax"].set_xlim(xmin, xmax)
    figure["ax"].set_ylim(ymin, ymax)

    return Plot(x=xvals,
                layers=to_render,
                fig=figure["fig"],
                ax=figure["ax"],
                filename=filename)
