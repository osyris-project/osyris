# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import importlib.metadata

from .config import config
from .core import Array, Datagroup, Dataset, Plot, Vector
from .plot import histogram1d, histogram2d, map, plot, scatter
from .spatial import extract_box, extract_sphere
from .units import units

try:
    __version__ = importlib.metadata.version(__package__ or __name__)
except importlib.metadata.PackageNotFoundError:
    __version__ = "0.0.0"

del importlib

__all__ = [
    "Array",
    "Datagroup",
    "Dataset",
    "Plot",
    "Vector",
    "config",
    "units",
    "histogram1d",
    "histogram2d",
    "scatter",
    "map",
    "plot",
    "extract_box",
    "extract_sphere",
]
