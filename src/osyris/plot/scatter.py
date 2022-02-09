# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

from pint.quantity import Quantity
from ..core import Plot, Array
from .. import units
from .render import render
from .parser import parse_layer


def scatter(x: Array,
            y: Array,
            logx: bool = False,
            logy: bool = False,
            loglog: bool = False,
            norm: str = None,
            filename: str = None,
            title: str = None,
            xmin: float = None,
            xmax: float = None,
            ymin: float = None,
            ymax: float = None,
            vmin: float = None,
            vmax: float = None,
            ax: object = None,
            **kwargs) -> Plot:
    """
    Make a 2D scatter plot with two variables as input.

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

    xvals = x.norm.values
    yvals = y.norm.values

    _, _, params = parse_layer(layer=None, norm=norm, vmin=vmin, vmax=vmax, **kwargs)
    to_render = [{
        "data": None,
        "mode": "scatter",
        "params": params,
    }]
    if "c" in params and not isinstance(params["c"], str):
        to_render[0].update({"unit": params["c"].unit.units, "name": params["c"].name})
        params["c"] = params["c"].norm.values
    if "s" in params:
        unit = None
        if isinstance(params["s"], Array):
            unit = params["s"].unit.units
        if isinstance(params["s"], Quantity):
            unit = params["s"].units
        if unit is not None:
            if unit != units.dimensionless:
                if x.unit.units != y.unit.units:
                    raise RuntimeError("Scatter: an Array with a unit was supplied "
                                       "as a size, but the units of the x and y "
                                       "Arrays do not agree. The size must either "
                                       "be a float or a dimensionless Array.")
                params["s"] = params["s"].to(x.unit.units)

    figure = render(x=xvals, y=yvals, data=to_render, logx=logx, logy=logy, ax=ax)

    figure["ax"].set_xlabel(x.label)
    figure["ax"].set_ylabel(y.label)
    figure["ax"].set_xlim(xmin, xmax)
    figure["ax"].set_ylim(ymin, ymax)

    return Plot(x=xvals,
                y=yvals,
                layers=to_render,
                fig=figure["fig"],
                ax=figure["ax"],
                filename=filename)
