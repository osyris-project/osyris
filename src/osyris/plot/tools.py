# # SPDX-License-Identifier: BSD-3-Clause
# # Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)

# from matplotlib.colors import LogNorm, Normalize, SymLogNorm

# from ..core.layer import Layer


# # def set_layer_norm(layer: Layer):
# #     print("layer.norm", layer.norm)
# #     if layer.norm is None:
# #         layer.norm = Normalize(vmin=layer.vmin, vmax=layer.vmax)
# #     if isinstance(layer.norm, str):
# #         norm_lowercase = layer.norm.lower()
# #         if norm_lowercase == "log":
# #             layer.norm = LogNorm(vmin=layer.vmin, vmax=layer.vmax)
# #         elif norm_lowercase == "symlog":
# #             layer.norm = SymLogNorm(
# #                 linthresh=1e-2, vmin=layer.vmin, vmax=layer.vmax, base=10
# #             )
# #         elif norm_lowercase == "linear":
# #             layer.norm = Normalize(vmin=layer.vmin, vmax=layer.vmax)
# #         else:
# #             raise RuntimeError(
# #                 f"Unknown norm keyword '{layer.norm}'.\nAvailable keywords"
# #                 " are 'log', 'symlog' and 'linear'."
# #             )

# #     layer.kwargs.update(norm=layer.norm)


# def parse_layer(
#     layer: Layer, mode=None, operation=None, norm=None, vmin=None, vmax=None, **kwargs
# ):
#     out = layer.copy()
#     if out.mode is None:
#         out.mode = mode
#     if out.operation is None:
#         out.operation = operation
#     if out.norm is None:
#         out.norm = norm
#     if out.vmin is None:
#         out.vmin = vmin
#     if out.vmax is None:
#         out.vmax = vmax
#     out.kwargs.update(
#         {key: value for key, value in kwargs.items() if key not in out.kwargs}
#     )

#     # Set the norm
#     if out.norm is None:
#         out.norm = Normalize(vmin=out.vmin, vmax=out.vmax)
#     if isinstance(out.norm, str):
#         norm_lowercase = out.norm.lower()
#         if norm_lowercase == "log":
#             out.norm = LogNorm(vmin=out.vmin, vmax=out.vmax)
#         elif norm_lowercase == "symlog":
#             out.norm = SymLogNorm(linthresh=1e-2, vmin=out.vmin, vmax=out.vmax, base=10)
#         elif norm_lowercase == "linear":
#             out.norm = Normalize(vmin=out.vmin, vmax=out.vmax)
#         else:
#             raise RuntimeError(
#                 f"Unknown norm keyword '{out.norm}'.\nAvailable keywords"
#                 " are 'log', 'symlog' and 'linear'."
#             )
#     out.kwargs.update(norm=out.norm)
#     return out
