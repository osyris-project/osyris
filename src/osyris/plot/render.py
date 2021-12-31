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
        "image": "pcolormesh",
        "imshow": "pcolormesh"
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

        if func == "line_integral_convolution" and "color" in item["params"]:
            cblabel = make_label(name=item["params"]["color"].name,
                                 unit=item["params"]["color"].unit.units)
        else:
            cblabel = make_label(name=item.get("name", ""), unit=item.get("unit", ""))

        mpl_objects.append(
            getattr(wrappers, func)(ax=ax,
                                    x=x,
                                    y=y,
                                    z=item["data"],
                                    cbar=cbar,
                                    cblabel=cblabel,
                                    **item["params"]))

    out["objects"] = mpl_objects
    return out
