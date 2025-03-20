# SPDX-License-Identifier: BSD-3-Clause

"""Osyris: A Python package for the analysis of astrophysical simulations

isort:skip_file
"""

import importlib.metadata

from .config import config
from .units import units
from .core import Array, Datagroup, Dataset, Plot, Vector, VectorBasis
from .io import RamsesDataset
from .plot import hist1d, hist2d, map, plot, scatter
from .spatial import extract_box, extract_sphere

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
    "VectorBasis",
    "RamsesDataset",
    "config",
    "units",
    "hist1d",
    "hist2d",
    "scatter",
    "map",
    "plot",
    "extract_box",
    "extract_sphere",
]
