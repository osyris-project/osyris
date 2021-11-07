# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

from matplotlib.colors import LogNorm, Normalize


def get_norm(norm=None, vmin=None, vmax=None):
    if norm == "log":
        return LogNorm(vmin=vmin, vmax=vmax)
    else:
        return Normalize(vmin=vmin, vmax=vmax)


def parse_layer(layer,
                mode=None,
                norm=None,
                vmin=None,
                vmax=None,
                operation=None,
                **kwargs):

    if isinstance(layer, dict):
        params = {
            key: layer[key]
            for key in set(layer.keys()) - set(["data", "mode", "operation"])
        }
        if "norm" not in params:
            params["norm"] = norm
        if "vmin" in params:
            vmin = params["vmin"]
            del params["vmin"]
        if "vmax" in params:
            vmax = params["vmax"]
            del params["vmax"]

        params["norm"] = get_norm(norm=params["norm"], vmin=vmin, vmax=vmax)

        for key, arg in kwargs.items():
            if key not in params:
                params[key] = arg

        settings = {}
        for key in ["mode", "operation"]:
            settings[key] = layer[key] if key in layer else eval(key)
        return layer["data"], settings, params
    else:
        params = {"norm": get_norm(norm=norm, vmin=vmin, vmax=vmax)}
        settings = {"mode": mode, "operation": operation}
        params.update(kwargs)
        return layer, settings, params
