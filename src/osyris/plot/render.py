# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

import matplotlib.pyplot as plt

from . import wrappers
from .. import config
from ..core.tools import make_label


def render(x=None, y=None, data=None, logx=False, logy=False, ax=None):
    """
    Use matplotlib to plot histogram, slice or column density maps
    """

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    out = {"fig": fig, "ax": ax}
    if data is None:
        return out

    function_map = {
        "vec": "quiver",
        "vector": "quiver",
        "stream": "streamplot",
        "lic": "line_integral_convolution",
        None: config.parameters["render_mode"],
        "image": "contourf",
        "imshow": "contourf"
    }

    mpl_objects = []
    for item in data:
        func = item["mode"]
        if func in function_map:
            func = function_map[func]

        if "cbar" in item["params"]:
            cbar = item["params"]["cbar"]
            del item["params"]["cbar"]
        else:
            cbar = True

        mpl_objects.extend(
            getattr(wrappers, func)(ax, x, y, item["data"], **item["params"]))

        need_cbar = False

        ind_render = -1
        name = item["name"]
        unit = item["unit"]

        if func == "line_integral_convolution" and "color" in item["params"]:
            need_cbar = True
            ind_render = -2
            name = item["params"]["color"].name
            unit = item["params"]["color"].unit.units

        if func in ["contourf", "pcolormesh"]:
            need_cbar = True
        if (func == "scatter") and ("c" in item["params"]):
            if not isinstance(item["params"]["c"], str):
                need_cbar = True
                name = item["name"]
                unit = item["unit"]
        if need_cbar and cbar:
            cb = plt.colorbar(mpl_objects[ind_render], ax=ax, cax=None)
            cb.set_label(make_label(name=name, unit=unit))
            cb.ax.yaxis.set_label_coords(-1.1, 0.5)
    out["objects"] = mpl_objects
    return out
