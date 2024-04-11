# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)

from matplotlib.colors import LogNorm, Normalize, SymLogNorm

from ..core.layer import Layer


def set_layer_norm(layer: Layer):
    print("layer.norm", layer.norm)
    if layer.norm is None:
        layer.norm = Normalize(vmin=layer.vmin, vmax=layer.vmax)
    if isinstance(layer.norm, str):
        norm_lowercase = layer.norm.lower()
        if norm_lowercase == "log":
            layer.norm = LogNorm(vmin=layer.vmin, vmax=layer.vmax)
        elif norm_lowercase == "symlog":
            layer.norm = SymLogNorm(
                linthresh=1e-2, vmin=layer.vmin, vmax=layer.vmax, base=10
            )
        elif norm_lowercase == "linear":
            layer.norm = Normalize(vmin=layer.vmin, vmax=layer.vmax)
        else:
            raise RuntimeError(
                f"Unknown norm keyword '{layer.norm}'.\nAvailable keywords"
                " are 'log', 'symlog' and 'linear'."
            )

    layer.kwargs.update(norm=layer.norm)
