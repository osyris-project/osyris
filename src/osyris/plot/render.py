# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import matplotlib.pyplot as plt
from ..config import parameters
from . import wrappers
from ..core.tools import make_label


def render(x, y, data, logx=False, logy=False, ax=None):
    """
    Use matplotlib to plot histogram, slice or column density maps
    """

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    function_map = {
        "vec": "quiver",
        "vector": "quiver",
        "stream": "streamplot",
        None: parameters["render_mode"],
        "image": "pcolormesh",
        "imshow": "pcolormesh"
    }

    mpl_objects = []
    cbar = None
    for item in data:
        func = item["mode"]
        if func in function_map:
            func = function_map[func]

        mpl_objects.append(
            getattr(wrappers, func)(ax, x, y, item["data"], **item["params"]))

        if func in ["contourf", "pcolormesh"] and cbar is None:
            cbar = plt.colorbar(mpl_objects[-1], ax=ax, cax=None)
            cbar.set_label(make_label(name=item["name"], unit=item["unit"]))
            cbar.ax.yaxis.set_label_coords(-1.1, 0.5)

    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    return {"fig": fig, "ax": ax, "objects": mpl_objects}
