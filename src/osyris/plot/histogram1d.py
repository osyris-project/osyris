# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..core import Plot
from .render import render
from ..core.tools import to_bin_centers, finmin, finmax


def histogram1d(*layers,
                bins=50,
                weights=None,
                ax=None,
                logx=False,
                logy=False,
                loglog=False,
                filename=None,
                title=None,
                ymin=None,
                ymax=None,
                **kwargs):
    """
    Plot a 1D histogram with arbitrary number of variables as input.
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
                xedges = np.logspace(np.log10(xmin), np.log10(xmax), bins + 1)
            else:
                xedges = np.linspace(xmin, xmax, bins + 1)
        else:
            xedges = params['bins']

        ydata, _, _ = figure["ax"].hist(xvals,
                                        bins=xedges,
                                        weights=params['weights'],
                                        **extra_args)

        figure["ax"].set_xlabel(params['data'].label)

    figure["ax"].set_ylim(ymin, ymax)
    return Plot(x=to_bin_centers(xedges), y=ydata, fig=figure["fig"], ax=figure["ax"])
