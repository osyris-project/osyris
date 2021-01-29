
def quiver(ax, x, y, z, density=1, **kwargs):

    default_args = {
        "angles": "xy",
        "pivot": "mid",
        "color": "w"
        }

    default_args.update(kwargs)

    skip = int(round(4.0 / density))
    skip = (slice(None,None,skip),slice(None,None,skip))

    return ax.quiver(x[skip[0]], y[skip[1]], z[..., 0][skip], z[..., 1][skip],
        **default_args)


def pcolormesh(ax, x, y, z, **kwargs):

    default_args = {
        "shading": "nearest",
        }

    default_args.update(kwargs)

    return ax.pcolormesh(x, y, z, **default_args)


def contour(ax, x, y, z, labels=True, **kwargs):

    cs = ax.contour(x, y, z, **kwargs)
    if labels:
        ax.clabel(cs, inline=1, fontsize=10)
    return cs

def contourf(ax, x, y, z, **kwargs):

    return ax.contourf(x, y, z, **kwargs)


def streamplot(ax, x, y, z, **kwargs):

    default_args = {
        "color": "w"
        }

    default_args.update(kwargs)

    return ax.streamplot(x, y, z[..., 0], z[..., 1], **default_args)
