# SPDX-License-Identifier: BSD-3-Clause
from functools import partial

import numpy as np

from .tools import make_label


class Base:
    def __repr__(self):
        return str(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    @property
    def label(self):
        """
        Return a label for the object with name and unit.
        """
        return make_label(name=self.name, unit=self.unit)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Numpy array_ufunc protocol to allow Array to work with numpy ufuncs.
        """
        if method != "__call__":
            # Only handle ufuncs as callables
            return NotImplemented
        return self._wrap_numpy(ufunc, *inputs, **kwargs)

    def __array_function__(self, func, types, args, kwargs):
        """
        Numpy array_function protocol to allow Array to work with numpy
        functions.
        """
        return self._wrap_numpy(func, *args, **kwargs)

    def __getattr__(self, name):
        if hasattr(np, name):
            return partial(getattr(np, name), self)
        raise AttributeError(f"No attribute named '{name}'")
