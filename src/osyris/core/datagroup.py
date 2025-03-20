# SPDX-License-Identifier: BSD-3-Clause
import numpy as np

from .layer import Layer
from .tools import bytes_to_human_readable
from .vector import Vector


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
        """
        Create a shallow copy of the Datagroup.
        """
        return self.__class__(**{key: array for key, array in self.items()})

    def keys(self):
        """
        The keys of the Datagroup (iterable).
        """
        return self._container.keys()

    def items(self):
        """
        The items of the Datagroup (iterable).
        """
        return self._container.items()

    def values(self):
        """
        The values of the Datagroup (iterable).
        """
        return self._container.values()

    def nbytes(self):
        """
        The number of bytes used by the Datagroup.
        """
        return np.sum([item.nbytes for item in self.values()])

    def print_size(self):
        """
        Return the size of the Datagroup in human readable format.
        """
        return bytes_to_human_readable(self.nbytes())

    @property
    def shape(self):
        """
        The shape of the Datagroup.
        """
        if len(self) == 0:
            return ()
        else:
            return self[list(self.keys())[0]].shape

    def sortby(self, key):
        """
        Sort the Datagroup by key.

        Parameters
        ----------
        key : str or list
            The key to sort the Datagroup by. If a list of indices is given, the
            Datagroup will be sorted by the indices.
        """
        if key is not None:
            if isinstance(key, str):
                key = np.argsort(self[key]).values
            for var in self.keys():
                self[var] = self[var][key]

    def clear(self):
        """
        Clear the Datagroup.
        """
        self._container.clear()

    def get(self, key, default):
        """
        Get the value of a key in the Datagroup.
        """
        return self._container.get(key, default)

    def pop(self, key):
        """
        Pop a key from the Datagroup.
        """
        return self._container.pop(key)

    def update(self, *args, **kwargs):
        """
        Update the Datagroup with new values.
        """
        d = dict(*args, **kwargs)
        for key, value in d.items():
            self[key] = value

    def layer(self, key: str, **kwargs) -> Layer:
        """
        Make a layer for map plots which contains mesh information.
        """
        keys = ("position", "dx", "mass", "velocity")
        return Layer(
            data=self[key],
            aux={k: self[k] for k in keys if k in self},
            **kwargs,
        )

    def to_pandas(self):
        """
        Convert the Datagroup to a pandas DataFrame.
        """
        import pandas as pd

        data = {}
        for key, item in self.items():
            data[key] = item.norm.values
            if isinstance(item, Vector):
                for c, xyz in item._xyz.items():
                    data[f"{key}_{c}"] = xyz.values
        return pd.DataFrame(data)
