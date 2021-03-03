import numpy as np
from ..config import parameters


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
        kwargs["cmap"] = parameters["cmap"]
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
        kwargs["cmap"] = parameters["cmap"]
    return ax.contourf(x, y, z, **kwargs)


def streamplot(ax, x, y, z, **kwargs):
    """
    Wrapper around Matplotlib's streamplot plot.
    """
    default_args = {"color": "w"}
    default_args.update(kwargs)
    return ax.streamplot(x, y, z[..., 0], z[..., 1], **default_args)
