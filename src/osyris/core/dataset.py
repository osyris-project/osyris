# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
from .. import config
from ..io import Loader
from .datagroup import Datagroup
from .tools import bytes_to_human_readable


class Dataset:
    def __init__(self, nout=1, scale=None, path=""):
        self.groups = {}
        self.meta = {}

        if scale is None:
            scale = config.parameters["scale"]

        self.loader = Loader(nout=nout, scale=scale, path=path)
        self.meta.update(self.loader.load_metadata())

    def __iter__(self):
        return self.groups.__iter__()

    def __len__(self):
        return self.groups.__len__()

    def __getitem__(self, key):
        return self.groups[key]

    def __setitem__(self, key, value):
        if not isinstance(value, Datagroup):
            raise TypeError("Only objects of type Datagroup can be inserted into a "
                            "Dataset.")
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

    def keys(self):
        return self.groups.keys()

    def items(self):
        return self.groups.items()

    def values(self):
        return self.groups.values()

    def load(self, *args, **kwargs):
        groups = self.loader.load(*args, meta=self.meta, **kwargs)
        for name, group in groups.items():
            self[name] = group
        config.additional_variables(self)
        return self

    def nbytes(self):
        return np.sum([item.nbytes() for item in self.groups.values()])

    def print_size(self):
        return bytes_to_human_readable(self.nbytes())
