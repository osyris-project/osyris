# SPDX-License-Identifier: BSD-3-Clause

from matplotlib.colors import LogNorm, Normalize, SymLogNorm

from ..core.layer import Layer


def get_norm(norm=None, vmin=None, vmax=None):
    if norm is None:
        return Normalize(vmin=vmin, vmax=vmax)
    if isinstance(norm, str):
        norm_lowercase = norm.lower()
        if norm_lowercase == "log":
            return LogNorm(vmin=vmin, vmax=vmax)
        elif norm_lowercase == "symlog":
            return SymLogNorm(linthresh=1e-2, vmin=vmin, vmax=vmax, base=10)
        elif norm_lowercase == "linear":
            return Normalize(vmin=vmin, vmax=vmax)
        else:
            raise RuntimeError(
                "Unknown norm keyword '{}'.\nAvailable keywords"
                " are 'log', 'symlog' and 'linear'.".format(norm)
            )
    else:
        return norm


def parse_layer(
    layer: Layer,
    mode=None,
    operation=None,
    norm=None,
    vmin=None,
    vmax=None,
    bins=None,
    weights=None,
    **kwargs,
):
    out = layer.copy()
    if out.mode is None:
        out.mode = mode
    if out.operation is None:
        out.operation = operation
    if out.norm is None:
        out.norm = norm
    if out.vmin is None:
        out.vmin = vmin
    if out.vmax is None:
        out.vmax = vmax
    if out.bins is None:
        out.bins = bins
    if out.weights is None:
        out.weights = weights
    out.kwargs.update(
        {key: value for key, value in kwargs.items() if key not in out.kwargs}
    )
    return out
