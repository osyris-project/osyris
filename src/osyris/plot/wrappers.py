# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import numpy as np
from .. import config
from ..core import Array
from ..core.tools import to_bin_edges
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

    return [ax.quiver(*args, **default_args)]


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
    return [ax.pcolormesh(x, y, z, **default_args)]


def contour(ax, x, y, z, labels=True, **kwargs):
    """
    Wrapper around Matplotlib's contour plot.

    We add a small convenience argument `labels` that allows to simply add
    clabels with a simple `labels=True`.
    """
    cs = ax.contour(x, y, z, **kwargs)
    if labels:
        ax.clabel(cs, inline=1, fontsize=10)
    return [cs]


def contourf(ax, x, y, z, **kwargs):
    """
    Wrapper around Matplotlib's contourf plot.
    """
    if "cmap" not in kwargs:
        kwargs["cmap"] = config.parameters["cmap"]
    return [ax.contourf(x, y, z, **kwargs)]


def streamplot(ax, x, y, z, **kwargs):
    """
    Wrapper around Matplotlib's streamplot plot.
    """
    default_args = {"color": "w"}
    default_args.update(kwargs)
    return [ax.streamplot(x, y, z[..., 0], z[..., 1], **default_args)]


def scatter(ax, x, y, data, **kwargs):
    """
    Wrapper around Matplotlib's scatter plot.
    If a point size has a unit, use PatchCollection instead of scatter.
    We use the scatter API for the arguments.
    If PatchCollection is used, we convert the scatter args to the
    PatchCollection syntax (e.g. "c" -> "color").
    """
    default_args = {"c": "b", "edgecolors": "k"}
    default_args.update(kwargs)
    use_patchcollection = False
    if "s" in default_args:
        if isinstance(default_args["s"], Array):
            default_args["s"] = default_args["s"].norm.values
            use_patchcollection = True
        if isinstance(default_args["s"], Quantity):
            default_args["s"] = np.full_like(x, default_args["s"].magnitude)
            use_patchcollection = True
    if use_patchcollection:
        array = None
        norm = None
        if "c" in default_args:
            if "facecolors" not in default_args:
                if isinstance(default_args["c"], str):
                    default_args["facecolors"] = default_args["c"]
                else:
                    array = default_args["c"]
            del default_args["c"]
        if "norm" in default_args:
            norm = default_args["norm"]
            del default_args["norm"]
        try:
            iter(default_args["s"])
        except TypeError:
            default_args["s"] = [default_args["s"]] * len(x)
        patches = [
            plt.Circle([x_, y_], s) for x_, y_, s in zip(x, y, default_args["s"])
        ]
        del default_args["s"]
        coll = ax.add_collection(PatchCollection(patches, **default_args))
        if array is not None:
            coll.set_array(array)
        if norm is not None:
            coll.set_norm(norm)
        return [coll]
    else:
        return [ax.scatter(x, y, **default_args)]


def line_integral_convolution(ax, x, y, z, length=30, color=None, **kwargs):
    """
    Wrapper that plots a line integral convolution of a vector field.
    Uses alpha blending to merge the LIC with a user-requested color.
    """
    import lic as lic
    from matplotlib.cm import ScalarMappable

    xedges = to_bin_edges(x)
    yedges = to_bin_edges(y)

    # compute line integral convolution
    lic_res = lic.lic(z[..., 1], z[..., 0], length=length)

    plot_args = {**kwargs}
    plot_args["extent"] = [
        np.min(xedges), np.max(xedges),
        np.min(yedges), np.max(yedges)
    ]
    plot_args["origin"] = "lower"

    axes = []

    # plot color
    if color is not None:
        plot_args_color = {**plot_args}
        plot_args_color["alpha"] = 1.
        axes.append(ax.imshow(z[..., 2], **plot_args_color))

        plot_args["alpha"] = .3  # ready alpha blending

    # amplify contrast on lic
    lim = (.2, .5)
    lic_data_clip = np.clip(lic_res, lim[0], lim[1])
    lic_data_rgba = ScalarMappable(norm=None, cmap="binary").to_rgba(lic_data_clip)
    lic_data_clip_rescale = (lic_data_clip - lim[0]) / (lim[1] - lim[0])
    lic_data_rgba[..., 3] = lic_data_clip_rescale

    # plot the lic
    plot_args["cmap"] = "binary"
    axes.append(ax.imshow(lic_data_rgba, **plot_args))

    return axes
