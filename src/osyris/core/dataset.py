# SPDX-License-Identifier: BSD-3-Clause
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
        for key, group in dict(*args, **kwargs).items():
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
        header = "Dataset: "
        if "infile" in self.meta:
            header += f"{self.meta['infile']}: "
        header += f"{self.print_size()}\n\n"
        body = "\n\n".join([str(item) for item in self.values()])
        return header + body

    def __copy__(self):
        return self.copy()

    def copy(self):
        out = self.__class__(**dict(self.items()))
        out.meta = self.meta.copy()
        return out

    def keys(self):
        return self.groups.keys()

    def items(self):
        return self.groups.items()

    def values(self):
        return self.groups.values()

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

    def update(self, *args, **kwargs):
        d = dict(*args, **kwargs)
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
