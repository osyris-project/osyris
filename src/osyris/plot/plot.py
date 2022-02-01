# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from typing import Union
from ..core import Plot, Array
from .render import render


def plot(x: Array,
         *y: Union[Array, dict],
         logx: bool = False,
         logy: bool = False,
         loglog: bool = False,
         filename: str = None,
         title: str = None,
         xmin: float = None,
         xmax: float = None,
         ymin: float = None,
         ymax: float = None,
         legend=True,
         ax: object = None,
         **kwargs) -> Plot:
    """
    Make a 1D plot with two variables as input.

    This function has an API very close to that of matplotlib's ``scatter`` function.
    For the documentation of any parameters that are not listed below, see
    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.scatter.html.

    :param x: Array to use for x positions along the horizontal dimension.

    :param y: Array(s) to use for y positions along the vertical dimension.

    :param logx: If ``True``, use logarithmic scaling on the horizontal axis.
        Default is ``False``.

    :param logy: If ``True``, use logarithmic scaling on the vertical axis.
        Default is ``False``.

    :param loglog: If ``True``, use logarithmic scaling on the horizontal and
        vertical axes. Default is ``False``.

    :param filename: If specified, the returned figure is also saved to file.
        Default is ``None``.

    :param title: The title of the figure. Default is ``None``.

    :param xmin: Minimum value for the horizontal axis. Default is ``None``.

    :param xmax: Maximum value for the horizontal axis. Default is ``None``.

    :param ymin: Minimum value for the vertical axis. Default is ``None``.

    :param ymax: Maximum value for the vertical axis. Default is ``None``.

    :param ax: A matplotlib axes inside which the figure will be plotted.
        Default is ``None``, in which case some new axes a created.
    """
    if loglog:
        logx = logy = True

    to_plot = []
    yaxis_unit = None
    xvals = None

    if isinstance(x, dict):
        to_plot.append({
            "x": x["x"].norm,
            "y": x["y"].norm,
            "params": {
                **kwargs
            },
            "unit": x["y"].unit.units,
            "name": x["y"].name
        })
        yaxis_unit = x["y"].unit
    else:
        xvals = x.norm

    if isinstance(y, Array):
        y = [y]

    for layer in y:

        if xvals is not None:
            layer_x = xvals
        else:
            layer_x = to_plot[0]["x"]

        if isinstance(layer, dict):
            layer_dict = {
                "x": layer["x"].norm,
                "y": layer["y"].norm,
                "params": {
                    **kwargs
                },
                "unit": layer["y"].unit.units,
                "name": layer["y"].name
            }
            layer_dict["x"] = layer["x"].norm if "x" in layer else layer_x
            to_plot.append(layer_dict)
            layer_unit = layer["y"].unit
        else:
            to_plot.append({
                "x": layer_x,
                "y": layer.norm,
                "params": {
                    **kwargs
                },
                "unit": layer.unit.units,
                "name": layer.name
            })
            layer_unit = layer.unit

        if yaxis_unit is None:
            yaxis_unit = layer_unit
        else:
            if layer_unit != yaxis_unit:
                raise RuntimeError(
                    "Different layers in 1D plots must all have the same unit.")

    figure = render(logx=logx, logy=logy, ax=ax)
    for item in to_plot:
        sorting = np.argsort(item["x"].values)
        figure["ax"].plot(item["x"].values[sorting],
                          item['y'].values[sorting],
                          label=item["name"],
                          **item["params"])

    figure["ax"].set_xlabel(to_plot[0]['x'].label)
    figure["ax"].set_ylabel(to_plot[0]['y'].label)
    figure["ax"].set_xlim(xmin, xmax)
    figure["ax"].set_ylim(ymin, ymax)
    figure["ax"].set_title(title)
    if legend and len(to_plot) > 1:
        figure["ax"].legend()

    return Plot(layers=to_plot, fig=figure["fig"], ax=figure["ax"], filename=filename)
