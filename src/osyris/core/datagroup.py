# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np

from .layer import Layer
from .tools import bytes_to_human_readable


class Datagroup:
    def __init__(self, *args, **kwargs):
        self._container = {}
        self.name = ""
        for key, array in dict(*args, **kwargs).items():
            self[key] = array

    def __iter__(self):
        return self._container.__iter__()

    def __len__(self):
        return self._container.__len__()

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._container[key]
        else:
            d = self.__class__()
            for name, val in self.items():
                d[name] = val[key]
            return d

    def __setitem__(self, key, value):
        if self.shape and (self.shape != value.shape):
            raise ValueError(
                "Size mismatch on element insertion. Item "
                "shape is {} while container accepts shape {}.".format(
                    value.shape, self.shape
                )
            )
        value.name = key
        self._container[key] = value

    def __delitem__(self, key):
        return self._container.__delitem__(key)

    def __repr__(self):
        return str(self)

    def __str__(self):
        header = f"Datagroup: {self.name} {self.print_size()}\n"
        body = "\n".join([str(item) for item in self.values()])
        return header + body

    def __eq__(self, other):
        if self.keys() != other.keys():
            return False
        for key, value in self.items():
            if all(value != other[key]):
                return False
        return True

    def __copy__(self):
        return self.copy()

    def copy(self):
        return self.__class__(**{key: array for key, array in self.items()})

    def keys(self):
        return self._container.keys()

    def items(self):
        return self._container.items()

    def values(self):
        return self._container.values()

    def nbytes(self):
        return np.sum([item.nbytes for item in self.values()])

    def print_size(self):
        return bytes_to_human_readable(self.nbytes())

    @property
    def shape(self):
        if len(self) == 0:
            return ()
        else:
            return self[list(self.keys())[0]].shape

    def sortby(self, key):
        if key is not None:
            if isinstance(key, str):
                key = np.argsort(self[key]).values
            for var in self.keys():
                self[var] = self[var][key]

    def clear(self):
        self._container.clear()

    def get(self, key, default):
        return self._container.get(key, default)

    def pop(self, key):
        return self._container.pop(key)

    def update(self, *args, **kwargs):
        d = dict(*args, **kwargs)
        for key, value in d.items():
            self[key] = value

    def layer(self, key: str, **kwargs) -> Layer:
        """
        Make a layer for map plots which contains mesh information
        """
        keys = ("position", "dx", "mass", "velocity")
        return Layer(
            data=self[key],
            aux={k: self[k] for k in keys if k in self},
            **kwargs,
        )
