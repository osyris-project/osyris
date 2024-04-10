# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np

from .. import config
from ..units import UnitsLibrary, units
from .datagroup import Datagroup
from .tools import bytes_to_human_readable


class Dataset:
    def __init__(self, *args, **kwargs):
        self.groups = {}
        self.meta = {}
        self.loader = None
        self.units = None
        entries = kwargs
        if args:
            if len(args) > 1:
                raise TypeError("Only one positional argument is allowed.")
            entries.update(args[0])
        for key, group in entries.items():
            self[key] = group

    def __iter__(self):
        return self.groups.__iter__()

    def __len__(self):
        return self.groups.__len__()

    def __getitem__(self, key):
        return self.groups[key]

    def __setitem__(self, key, value):
        if not isinstance(value, Datagroup):
            raise TypeError(
                "Only objects of type Datagroup can be inserted into a " "Dataset."
            )
        self.groups.__setitem__(key, value)
        value.name = key
        value.parent = self

    def __delitem__(self, key):
        return self.groups.__delitem__(key)

    def __repr__(self):
        return str(self)

    def __str__(self):
        output = "Dataset: "
        if "infile" in self.meta:
            output += "{}: ".format(self.meta["infile"])
        output += "{}\n".format(self.print_size())
        for item in self.values():
            output += str(item) + "\n"
        return output

    def __copy__(self):
        return self.copy()

    def copy(self):
        out = self.__class__(
            **{key: group.copy() for key, group in self.groups.items()}
        )
        out.meta = self.meta.copy()
        return out

    def keys(self):
        return self.groups.keys()

    def items(self):
        return self.groups.items()

    def values(self):
        return self.groups.values()

    # def load(self, *args, **kwargs):
    #     groups = self.loader.load(*args, meta=self.meta, units=self.units, **kwargs)
    #     for name, group in groups.items():
    #         self[name] = group
    #     config.additional_variables(self)
    #     return self

    def nbytes(self):
        return np.sum([item.nbytes() for item in self.groups.values()])

    def print_size(self):
        return bytes_to_human_readable(self.nbytes())

    def clear(self):
        self.groups.clear()
        self.meta.clear()

    def get(self, key, default):
        return self.groups.get(key, default)

    def pop(self, key):
        return self.groups.pop(key)

    def update(self, d):
        for key, value in d.items():
            self[key] = value

    def set_units(self):
        self.units = UnitsLibrary(
            library=config.configure_units(
                units=units,
                unit_d=self.meta["unit_d"],
                unit_l=self.meta["unit_l"],
                unit_t=self.meta["unit_t"],
            ),
            default_unit=1.0 * units(""),
        )
