# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
import numpy as np
from .array import Array
from .tools import bytes_to_human_readable
# from ..utils import value_to_string


class Datagroup:
    def __init__(self):
        self._container = {}
        # self.meta = {}
        self.shape = None

    def __iter__(self):
        return self._container.__iter__()

    def __len__(self):
        return self._container.__len__()

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._container[key]
        else:
            d = Dataset()
            for name, val in self.items():
                d[name] = val[key]
            # d.meta.update(self.meta)
            return d

    def __setitem__(self, key, value):
        if len(value.shape) > 0:
            shape = value.shape[0]
        else:
            shape = 1
        if self.shape is None:
            self.shape = shape
        else:
            if self.shape != shape:
                raise RuntimeError(
                    "Size mismatch on element insertion. Item "
                    "shape is {} while container accepts shape {}.".format(
                        shape, self.shape))
        if isinstance(value, Array):
            value.name = key
            value.parent = self
            self._container[key] = value
        else:
            self._container[key] = Array(values=value, name=key, parent=self)

    def __delitem__(self, key):
        return self._container.__delitem__(key)

    def __repr__(self):
        return str(self)

    def __str__(self):
        output = "Dataset: "
        # if "infile" in self.meta:
        #     output += "{}: ".format(self.meta["infile"])
        output += "{}\n".format(self.print_size())
        for key, item in self.items():
            output += str(item) + "\n"
        return output

    def keys(self):
        return self._container.keys()

    def items(self):
        return self._container.items()

    def values(self):
        return self._container.values()

    def set_scale(self, scale):
        for key in ["x", "y", "z", "dx"]:
            if key in self:
                self[key].to(scale)

    def nbytes(self):
        return np.sum([item._array.nbytes for item in self.values()])

    def print_size(self):
        return bytes_to_human_readable(self.nbytes())
