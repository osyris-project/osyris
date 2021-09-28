# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

from pint.quantity import Quantity
from ..core import Plot, Array
from .. import units
from .render import render
from .parser import parse_layer


def scatter(x,
            y,
            ax=None,
            logx=False,
            logy=False,
            loglog=False,
            norm=None,
            filename=None,
            title=None,
            xmin=None,
            xmax=None,
            ymin=None,
            ymax=None,
            vmin=None,
            vmax=None,
            **kwargs):
    """
    Plot a 2D scatter plot with two variables as input.
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

    return Plot(x=xvals, y=yvals, layers=to_render, fig=figure["fig"], ax=figure["ax"])
