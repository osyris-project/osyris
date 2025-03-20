# SPDX-License-Identifier: BSD-3-Clause

from typing import Union

from pint import Quantity

from .. import units
from ..core import Array, Plot, Vector
from .parser import get_norm
from .render import render


def scatter(
    x: Array,
    y: Array,
    color: Union[str, Array] = None,
    size: Union[float, Array] = None,
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
    aspect: str = "auto",
    ax: object = None,
    **kwargs,
) -> Plot:
    """
    Make a 2D scatter plot with two variables as input.

    :param x: Array to use for scatter point positions along the horizontal dimension.

    :param y: Array to use for scatter point positions along the vertical dimension.

    :param color: The color of the scatter points. Can be a string or an Array. If an
        Array is supplied, a colormap is used. Default is ``None``.

    :param size: Size of the scatter points. If a float is provided, all points will
        have the same size. If an Array is provided, the size of each point will be
        determined by the values of the Array. Default is ``None``.

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

    :param aspect: The aspect ratio of the plot. Default is ``'auto'``.

    :param ax: A matplotlib axes inside which the figure will be plotted.
        Default is ``None``, in which case some new axes a created.
    """
    if loglog:
        logx = logy = True

    xvals = x.norm.values
    yvals = y.norm.values

    params = {k: v for k, v in kwargs.items() if k not in ["c", "s"]}
    to_render = {"data": None, "mode": "scatter"}

    if color is not None:
        if isinstance(color, str):
            params["c"] = color
        else:
            to_render.update({"unit": color.unit, "name": color.name})
            params.update(
                c=color.norm.values, norm=get_norm(norm=norm, vmin=vmin, vmax=vmax)
            )

    if size is not None:
        unit = None
        if isinstance(size, (Array, Vector)):
            unit = size.unit
        if isinstance(size, Quantity):
            unit = size.units
        if unit is not None:
            if unit != units("dimensionless"):
                if x.unit != y.unit:
                    raise RuntimeError(
                        "Scatter: an Array with a unit was supplied "
                        "as a size, but the units of the x and y "
                        "Arrays do not agree. The size must either "
                        "be a float or a dimensionless Array."
                    )
                size = size.to(x.unit)
        params["s"] = size

    to_render["params"] = params

    figure = render(x=xvals, y=yvals, data=[to_render], logx=logx, logy=logy, ax=ax)

    figure["ax"].set_xlabel(x.label)
    figure["ax"].set_ylabel(y.label)
    figure["ax"].set_title(title)
    figure["ax"].autoscale_view()
    _xmin, _xmax = figure["ax"].get_xlim()
    _ymin, _ymax = figure["ax"].get_ylim()
    if xmin is None:
        xmin = _xmin
    if xmax is None:
        xmax = _xmax
    if ymin is None:
        ymin = _ymin
    if ymax is None:
        ymax = _ymax
    figure["ax"].set_xlim(xmin, xmax)
    figure["ax"].set_ylim(ymin, ymax)
    figure["ax"].set_aspect(aspect)

    return Plot(
        x=xvals,
        y=yvals,
        layers=to_render,
        fig=figure["fig"],
        ax=figure["ax"],
        filename=filename,
    )
