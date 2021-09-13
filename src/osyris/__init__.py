# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

# flake8: noqa

from pint import UnitRegistry

units = UnitRegistry(system="cgs")

from .config import config
config.additional_units(units)

from .io import Loader
from .plot import histogram, plane
from .core import Array, Dataset, Plot
