# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from copy import copy, deepcopy
from .. import config
from ..io import Loader
from .datagroup import Datagroup
from .tools import bytes_to_human_readable
from ..units import units, UnitsLibrary


class Dataset:
    def __init__(self, nout=None, path=""):
        self.groups = {}
        self.meta = {}
        self.loader = None
        self.units = None
        if nout is not None:
            self.loader = Loader(nout=nout, path=path)
            self.meta.update(self.loader.load_metadata())
            self.set_units()
            self.meta["time"] *= self.units["time"]

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
        for key, item in self.items():
            output += str(item) + "\n"
        return output

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        nout = self.loader.nout if self.loader is not None else None
        path = self.loader.path if self.loader is not None else ""
        out = self.__class__(nout=nout, path=path)
        for key, group in self.groups.items():
            out[key] = deepcopy(group)
        out.meta = {key: copy(item) for key, item in self.meta.items()}
        return out

    def copy(self):
        nout = self.loader.nout if self.loader is not None else None
        path = self.loader.path if self.loader is not None else ""
        out = self.__class__(nout=nout, path=path)
        for key, group in self.groups.items():
            out[key] = group
        out.meta = self.meta.copy()
        return out

    def keys(self):
        return self.groups.keys()

    def items(self):
        return self.groups.items()

    def values(self):
        return self.groups.values()

    def load(self, *args, **kwargs):
        groups = self.loader.load(*args, meta=self.meta, units=self.units, **kwargs)
        for name, group in groups.items():
            self[name] = group
        config.additional_variables(self)
        return self

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
