# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)

from matplotlib.colors import LogNorm, Normalize, SymLogNorm

from ..core.layer import Layer


# def get_norm(norm=None, vmin=None, vmax=None):
#     if norm is None:
#         return Normalize(vmin=vmin, vmax=vmax)
#     if isinstance(norm, str):
#         norm_lowercase = norm.lower()
#         if norm_lowercase == "log":
#             return LogNorm(vmin=vmin, vmax=vmax)
#         elif norm_lowercase == "symlog":
#             return SymLogNorm(linthresh=1e-2, vmin=vmin, vmax=vmax, base=10)
#         elif norm_lowercase == "linear":
#             return Normalize(vmin=vmin, vmax=vmax)
#         else:
#             raise RuntimeError(
#                 "Unknown norm keyword '{}'.\nAvailable keywords"
#                 " are 'log', 'symlog' and 'linear'.".format(norm)
#             )
#     else:
#         return norm


# def parse_layer(
#     layer, mode=None, norm=None, vmin=None, vmax=None, operation=None, **kwargs
# ):
#     if isinstance(layer, dict):
#         params = {
#             key: layer[key]
#             for key in set(layer.keys()) - set(["data", "mode", "operation"])
#         }
#         if "norm" not in params:
#             params["norm"] = norm
#         if "vmin" in params:
#             vmin = params["vmin"]
#             del params["vmin"]
#         if "vmax" in params:
#             vmax = params["vmax"]
#             del params["vmax"]

#         params["norm"] = get_norm(norm=params["norm"], vmin=vmin, vmax=vmax)

#         for key, arg in kwargs.items():
#             if key not in params:
#                 params[key] = arg

#         settings = {}
#         for key in ["mode", "operation"]:
#             settings[key] = layer[key] if key in layer else eval(key)
#         return layer["data"], settings, params
#     else:
#         params = {"norm": get_norm(norm=norm, vmin=vmin, vmax=vmax)}
#         settings = {"mode": mode, "operation": operation}
#         params.update(kwargs)
#         return layer, settings, params


# def set_layer_norm(layer: Layer):
#     print("layer.norm", layer.norm)
#     if layer.norm is None:
#         layer.norm = Normalize(vmin=layer.vmin, vmax=layer.vmax)
#     if isinstance(layer.norm, str):
#         norm_lowercase = layer.norm.lower()
#         if norm_lowercase == "log":
#             layer.norm = LogNorm(vmin=layer.vmin, vmax=layer.vmax)
#         elif norm_lowercase == "symlog":
#             layer.norm = SymLogNorm(
#                 linthresh=1e-2, vmin=layer.vmin, vmax=layer.vmax, base=10
#             )
#         elif norm_lowercase == "linear":
#             layer.norm = Normalize(vmin=layer.vmin, vmax=layer.vmax)
#         else:
#             raise RuntimeError(
#                 f"Unknown norm keyword '{layer.norm}'.\nAvailable keywords"
#                 " are 'log', 'symlog' and 'linear'."
#             )

#     layer.kwargs.update(norm=layer.norm)


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

    # Set the norm
    if out.norm is None:
        out.norm = Normalize(vmin=out.vmin, vmax=out.vmax)
    if isinstance(out.norm, str):
        norm_lowercase = out.norm.lower()
        if norm_lowercase == "log":
            out.norm = LogNorm(vmin=out.vmin, vmax=out.vmax)
        elif norm_lowercase == "symlog":
            out.norm = SymLogNorm(linthresh=1e-2, vmin=out.vmin, vmax=out.vmax, base=10)
        elif norm_lowercase == "linear":
            out.norm = Normalize(vmin=out.vmin, vmax=out.vmax)
        else:
            raise RuntimeError(
                f"Unknown norm keyword '{out.norm}'. Available keywords"
                " are 'log', 'symlog' and 'linear'."
            )
    out.kwargs.update(norm=out.norm)
    return out
