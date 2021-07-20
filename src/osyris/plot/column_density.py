# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from .plane import plane
from ..core import Array, Plot
from .render import render


def column_density(*layers,
                   dz=None,
                   origin=None,
                   resolution=256,
                   ax=None,
                   plot=True,
                   **kwargs):
    """
    Plot a 2D slice through the data domain.
    """

    if isinstance(layers, Array):
        layers = [layers]

    dataset = layers[0].parent

    if not hasattr(dz, 'unit'):
        dz *= dataset["xyz"].unit
    # if dy is not None and not isinstance(dy, Quantity):
    #     dy *= dataset["xyz"].unit

    nsteps = 50

    step = dz / nsteps

    if origin is None:
        origin = Array(values=np.zeros([1, dataset.meta["ndim"]]),
                       unit=dataset["xyz"].unit)

    to_return = {}
    # to_render = [np.zeros([resolution] * 2)]
    to_render = None

    for i in range(nsteps):
        p = plane(*layers,
                  origin=origin - (0.5 * dz) + (step * i),
                  resolution=resolution,
                  ax=ax,
                  plot=False,
                  **kwargs)
        # to_return.update({"x": p.x, "y": p.y})
        if to_render is None:
            to_render = p.layers
        else:
            to_render[0]["data"] += p.layers[0]["data"]

    to_render[0]["data"] *= dz

    to_return = {"x": p.x, "y": p.y, "layers": to_render}
    if plot:
        # Render the map
        figure = render(x=p.x, y=p.y, data=to_render, ax=ax)
        figure["ax"].set_xlabel(dataset["xyz"].x.label)
        figure["ax"].set_ylabel(dataset["xyz"].y.label)
        if ax is None:
            figure["ax"].set_aspect("equal")
        to_return.update({"fig": figure["fig"], "ax": figure["ax"]})

    return Plot(**to_return)
