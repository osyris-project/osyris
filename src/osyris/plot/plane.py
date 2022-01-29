# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

from .map import map as _map
import warnings


def plane(*args, **kwargs):
    """
    Old alias for think map, will be deprecated soon.
    """
    warnings.warn("The plane function will be deprecated soon, use map instead.")
    return _map(*args, dz=None, **kwargs)
