import numpy as np
from matplotlib.colors import LogNorm, Normalize


def get_finite_inds(x):
    """
    Find indices of finite numbers.
    """
    return np.where(np.isfinite(x))[0]

def finmin(x):
    """
    Finite minimum.
    """
    return np.amin(x.take(get_finite_inds(x)))

def finmax(x):
    """
    Finite maximum.
    """
    return np.amax(x.take(get_finite_inds(x)))

def to_bin_centers(x):
    """
    Convert array edges to centers
    """
    return 0.5 * (x[1:] + x[:-1])

def to_bin_edges(x):
    """
    Convert array centers to edges
    """
    centers = to_bin_centers(x)
    left = centers[0] - (x[1] - x[0])
    right = centers[-1] + (x[-1] - x[-2])
    return np.concatenate(np.concatenate(left, center), right)

def get_norm(norm=None, vmin=None, vmax=None):
    if norm is not None:
        func = LogNorm if norm == "log" else Normalize
        norm = func(vmin=vmin, vmax=vmax)
    return norm

def parse_layer(entry, mode=None, norm=None, vmin=None, vmax=None, **kwargs):
    mode = "contourf" if mode is None else mode
    
    # print("NORM IS", norm)

    if isinstance(entry, dict):
        params = {key: entry[key] for key in set(entry.keys() - set(["data", "mode"]))}
        if "norm" not in params:
            params["norm"] = norm
        if "vmin" not in params:
            params["vmin"] = vmin
        if "vmax" not in params:
            params["vmax"] = vmax

        params["norm"] = get_norm(norm=params["norm"],
            vmin=params["vmin"], vmax=params["vmax"])

        for key, arg in kwargs.items():
            if key not in params:
                params[key] = arg
        # print(params)
        # params = {key: entry[key] for key in set(entry.keys() - set(["data"]))}

        if "mode" in entry:
            layer_mode = entry["mode"]
        else:
            layer_mode = mode

        return entry["data"], layer_mode, params
        # params = {key: entry[key] for key in set(entry.keys() - set(["data"]))}
        # return entry["data"], kwargs
        # ret
    else:
        params = {"mode": mode,
                   "norm": get_norm(norm=norm, vmin=vmin, vmax=vmax),
                   "vmin": vmin,
                   "vmax": vmax}
        params.update(kwargs)
        return entry, params
