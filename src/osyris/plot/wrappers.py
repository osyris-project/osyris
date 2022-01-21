# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
from .. import config
from ..core import Array
from contextlib import redirect_stderr
import io
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
import numpy as np
from pint.quantity import Quantity


def _add_colorbar(obj, ax, cax=None, label=None):
    """
    Add colorbar to plot for a given artist.
    """
    cb = plt.colorbar(obj, ax=ax, cax=cax)
    cb.set_label(label)
    cb.ax.yaxis.set_label_coords(-1.05, 0.5)


def quiver(ax,
           x,
           y,
           z,
           cbar=False,
           cblabel=None,
           density=1,
           color="w",
           zorder=2,
           **kwargs):
    """
    Wrapper around Matplotlib's quiver plot.

    We have a special argument `density` which aims to control the density of
    arrows in a similar manner to the density of streamlines in the streamplot.
    Matplotlib's own way of controlling the density of arrows is not so
    obvious.
    """
    default_args = {"angles": "xy", "pivot": "mid", "zorder": zorder}
    default_args.update(kwargs)

    skips = np.around(np.array(z.shape) * 4.0 / 128.0 / density).astype(np.int)
    skip = (slice(None, None, skips[0]), slice(None, None, skips[1]))

    args = [x[skip[0]], y[skip[1]], z[..., 0][skip], z[..., 1][skip]]
    if isinstance(color, str):
        default_args["color"] = color
    else:
        args.append(z[..., 2][skip])

    out = ax.quiver(*args, **default_args)
    if cbar and not isinstance(color, str):
        _add_colorbar(obj=out, ax=ax, label=cblabel)
    return out


def pcolormesh(ax, x, y, z, cbar=False, cblabel=None, zorder=1, **kwargs):
    """
    Wrapper around Matplotlib's pcolormesh plot.
    """
    default_args = {
        "shading": "nearest",
        "zorder": zorder,
        "cmap": config.parameters["cmap"]
    }
    default_args.update(kwargs)
    out = ax.pcolormesh(x, y, z, **default_args)
    if cbar:
        _add_colorbar(obj=out, ax=ax, label=cblabel)
    return out


def contour(ax,
            x,
            y,
            z,
            cbar=False,
            cblabel=None,
            labels=True,
            zorder=2,
            fmt='%1.3f',
            **kwargs):
    """
    Wrapper around Matplotlib's contour plot.

    We add a small convenience argument `labels` that allows to simply add
    clabels with a simple `labels=True`.
    """
    cs = ax.contour(x, y, z, zorder=zorder, **kwargs)
    if labels:
        ax.clabel(cs, cs.levels, inline=1, fontsize=10, fmt=fmt)
    return cs


def contourf(ax, x, y, z, cbar=False, cblabel=None, zorder=1, **kwargs):
    """
    Wrapper around Matplotlib's contourf plot.
    """
    if "cmap" not in kwargs:
        kwargs["cmap"] = config.parameters["cmap"]
    out = ax.contourf(x, y, z, **kwargs)
    if cbar:
        _add_colorbar(obj=out, ax=ax, label=cblabel)
    return out


def streamplot(ax, x, y, z, cbar=False, cblabel=None, color='w', zorder=2, **kwargs):
    """
    Wrapper around Matplotlib's streamplot plot.
    """
    default_args = {"color": "w", "zorder": zorder, "cmap": config.parameters["cmap"]}
    default_args.update(kwargs)
    if isinstance(color, str):
        default_args["color"] = color
    else:
        default_args["color"] = z[..., 2]
        if default_args["norm"].vmin is None:
            default_args["norm"].vmin = default_args["color"].min()
        if default_args["norm"].vmax is None:
            default_args["norm"].vmax = default_args["color"].max()
    out = ax.streamplot(x, y, z[..., 0], z[..., 1], **default_args)
    if cbar and not isinstance(color, str):
        _add_colorbar(obj=out.lines, ax=ax, label=cblabel)
    return out


def scatter(ax, x, y, z, cbar=False, cblabel=None, zorder=2, **kwargs):
    """
    Wrapper around Matplotlib's scatter plot.
    If a point size has a unit, use PatchCollection instead of scatter.
    We use the scatter API for the arguments.
    If PatchCollection is used, we convert the scatter args to the
    PatchCollection syntax (e.g. "c" -> "color").
    """
    default_args = {"c": "b", "edgecolors": "k", "zorder": zorder}
    default_args.update(kwargs)
    use_patchcollection = False
    need_cbar = False
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
            need_cbar = not isinstance(default_args["c"], str)
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
        out = ax.add_collection(PatchCollection(patches, **default_args))
        if array is not None:
            out.set_array(array)
        if norm is not None:
            out.set_norm(norm)
    else:
        need_cbar = not isinstance(default_args.get("c", ""), str)
        out = ax.scatter(x, y, **default_args)
    if cbar and need_cbar:
        _add_colorbar(obj=out, ax=ax, label=cblabel)
    return out


def line_integral_convolution(ax,
                              x,
                              y,
                              z,
                              cbar=False,
                              cblabel=None,
                              length=None,
                              color=None,
                              **kwargs):
    """
    Wrapper that plots a line integral convolution of a vector field.
    Uses alpha blending to merge the LIC with a user-requested color.
    """
    from lic import lic

    # Compute line integral convolution
    if length is None:
        length = int(max(z.shape[:-1]) * 15 / 128)
    with redirect_stderr(io.StringIO()) as _:
        lic_res = lic(z[..., 1], z[..., 0], length=length)

    if color is not None:
        plot_args = {**kwargs}
        base_args = {"norm": None, "cmap": None}
        for key in base_args:
            if key in plot_args:
                base_args[key] = plot_args[key]
                del plot_args[key]
        scalar_map = ScalarMappable(**base_args)
        base_rgba = scalar_map.to_rgba(z[..., 2])
        base_alpha = 1.0

        lim = (.2, .5)
        lic_data_clip = np.clip(lic_res, lim[0], lim[1])
        lic_rgba = ScalarMappable(norm=None, cmap="binary").to_rgba(lic_data_clip)
        lic_data_clip_rescale = (lic_data_clip - lim[0]) / (lim[1] - lim[0])
        lic_alpha = lic_data_clip_rescale.reshape(lic_data_clip_rescale.shape +
                                                  (1, )) * 0.3
        # Perform alpha blending manually
        rgba = lic_rgba * lic_alpha + (base_alpha * (1.0 - lic_alpha)) * base_rgba
        out = ax.imshow(rgba,
                        extent=[
                            0.5 * (3 * x[0] - x[1]), 0.5 * (3 * x[-1] - x[-2]),
                            0.5 * (3 * y[0] - y[1]), 0.5 * (3 * y[-1] - y[-2])
                        ],
                        origin='lower',
                        zorder=1)
        # Add the colorbar using the ScalarMappable
        scalar_map.set_array(z[..., 2])
        if cbar:
            _add_colorbar(obj=scalar_map, ax=ax, label=cblabel)
        return out
    else:
        return pcolormesh(ax=ax, x=x, y=y, z=lic_res, cbar=cbar, cblabel="", **kwargs)
