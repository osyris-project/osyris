# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from typing import Union, Iterable
from ..core import Plot, Array
from .render import render
from ..core.tools import to_bin_centers, finmin, finmax


def histogram1d(*layers: Union[Iterable, Array],
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
                **kwargs) -> Plot:
    """
    Plot a 1D histogram with arbitrary number of variables as input.
    When a vector quantity is supplied, the function will histogram the norm of
    the vectors.


    This function has an API very close to that of matplotlib's ``hist`` function.
    For the documentation of any parameters that are not listed below, see
    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html.

    :param layers: Dicts or Arrays representing the quantities to be mapped onto the
        colormap of the generated image.

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
        if isinstance(layer, dict):
            params = {}
            extra_args = {}
            for key, param in layer.items():
                if key in ["data", "bins", "weights"]:
                    params[key] = param
                else:
                    extra_args[key] = param
            for key, arg in {'bins': bins, 'weights': weights}.items():
                if key not in params:
                    params[key] = arg
        else:
            params = {'data': layer, 'bins': bins, 'weights': weights}
            extra_args = kwargs

        xvals = params['data'].norm.values
        if params['weights'] is not None:
            params['weights'] = params['weights'].norm.values

        # Construct some bin edges
        if isinstance(params['bins'], int):
            xmin = finmin(xvals)
            xmax = finmax(xvals)
            if logx:
                xedges = np.logspace(np.log10(xmin), np.log10(xmax), params['bins'] + 1)
            else:
                xedges = np.linspace(xmin, xmax, params['bins'] + 1)
        else:
            xedges = params['bins']

        ydata, _, _ = figure["ax"].hist(xvals,
                                        bins=xedges,
                                        weights=params['weights'],
                                        **extra_args)

        figure["ax"].set_xlabel(params['data'].label)

    figure["ax"].set_ylim(ymin, ymax)
    return Plot(x=to_bin_centers(xedges),
                y=ydata,
                fig=figure["fig"],
                ax=figure["ax"],
                filename=filename)
