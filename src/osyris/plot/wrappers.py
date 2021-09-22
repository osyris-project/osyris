# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import numpy as np
from .. import config
from ..core import Array
from pint.quantity import Quantity


def quiver(ax, x, y, z, density=1, color="w", **kwargs):
    """
    Wrapper around Matplotlib's quiver plot.

    We have a special argument `density` which aims to control the density of
    arrows in a similar manner to the density of streamlines in the streamplot.
    Matplotlib's own way of controlling the density of arrows is not so
    obvious.
    """
    default_args = {
        "angles": "xy",
        "pivot": "mid",
    }
    default_args.update(kwargs)

    skips = np.around(np.array(z.shape) * 4.0 / 128.0 / density).astype(np.int)
    skip = (slice(None, None, skips[0]), slice(None, None, skips[1]))

    args = [x[skip[0]], y[skip[1]], z[..., 0][skip], z[..., 1][skip]]
    if isinstance(color, str):
        default_args["color"] = color
    else:
        args.append(z[..., 2][skip])

    return ax.quiver(*args, **default_args)


def pcolormesh(ax, x, y, z, **kwargs):
    """
    Wrapper around Matplotlib's pcolormesh plot.
    """
    default_args = {
        "shading": "nearest",
    }
    default_args.update(kwargs)
    if "cmap" not in kwargs:
        kwargs["cmap"] = config.parameters["cmap"]
    return ax.pcolormesh(x, y, z, **default_args)


def contour(ax, x, y, z, labels=True, **kwargs):
    """
    Wrapper around Matplotlib's contour plot.

    We add a small convenience argument `labels` that allows to simply add
    clabels with a simple `labels=True`.
    """
    cs = ax.contour(x, y, z, **kwargs)
    if labels:
        ax.clabel(cs, inline=1, fontsize=10)
    return cs


def contourf(ax, x, y, z, **kwargs):
    """
    Wrapper around Matplotlib's contourf plot.
    """
    if "cmap" not in kwargs:
        kwargs["cmap"] = config.parameters["cmap"]
    return ax.contourf(x, y, z, **kwargs)


def streamplot(ax, x, y, z, **kwargs):
    """
    Wrapper around Matplotlib's streamplot plot.
    """
    default_args = {"color": "w"}
    default_args.update(kwargs)
    return ax.streamplot(x, y, z[..., 0], z[..., 1], **default_args)


def scatter(ax, x, y, data, **kwargs):
    """
    Wrapper around Matplotlib's scatter plot.
    If a point size has a unit, use PatchCollection instead of scatter.
    """
    default_args = {"c": "b"}
    default_args.update(kwargs)
    use_patchcollection = False
    if "s" in kwargs:
        if isinstance(kwargs["s"], Array):
            kwargs["s"] = kwargs["s"].norm.values
            use_patchcollection = True
        if isinstance(kwargs["s"], Quantity):
            kwargs["s"] = kwargs["s"].magnitude
            use_patchcollection = True
    if use_patchcollection:
        patches = [plt.Circle([x_, y_], s) for x_, y_, s in zip(x, y, kwargs["s"])]
        del kwargs["s"]
        array = None
        if "c" in kwargs:
            kwargs["array"] = kwargs["c"]
            del kwargs["c"]
        return ax.add_collection(PatchCollection(patches, **kwargs))
    else:
        return ax.scatter(x, y, **default_args)
