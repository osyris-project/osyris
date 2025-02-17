# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)

from pint import Quantity, Unit, UnitRegistry
from pint import compat

from .. import config


class Units:
    def __init__(self):
        self._ureg = UnitRegistry(system="cgs")
        config.configure_constants(self._ureg)

    def __call__(self, arg):
        if isinstance(arg, Quantity):
            raise TypeError("Cannot create unit from a Quantity. Use `.units` instead.")
        if isinstance(arg, Unit):
            return arg
        return self._ureg(arg).units

    def define(self, *args, **kwargs):
        self._ureg.define(*args, **kwargs)


units = Units()

# Add compatibility for upcasting
compat.upcast_type_map["osyris.core.array.Array"] = None
compat.upcast_type_map["osyris.core.array.Vector"] = None
