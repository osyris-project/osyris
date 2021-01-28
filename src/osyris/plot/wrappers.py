
def quiver(ax, x, y, u, v, density=1, **kwargs):

    default_args = {
        "angles": "xy",
        "pivot": "mid",
        "color": "w"
        }

    default_args.update(kwargs)

    skip = int(round(4.0 / density))
    skip = (slice(None,None,skip),slice(None,None,skip))

    return ax.quiver(x[skip[0]], y[skip[1]], u[skip], v[skip],
        **default_args)


def pcolormesh(ax, x, y, z, **kwargs):

    return ax.pcolormesh(x, y, z, **kwargs)


def contour(ax, x, y, z, **kwargs):

    return ax.contour(x, y, z, **kwargs)

def contourf(ax, x, y, z, **kwargs):

    return ax.contourf(x, y, z, **kwargs)


def streamplot(ax, x, y, u, v, **kwargs):

    default_args = {
        "color": "w"
        }

    default_args.update(kwargs)

    return ax.streamplot(x, y, u, v, **default_args)
