# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
from .array import Array
from .operators import add, sub
from .tools import value_to_string, make_label
from .. import units


class Vector(Array):
    def __getitem__(self, slice_):
        slice_ = tuple((slice_, )) + (slice(None, None, None), )
        return self.__class__(values=self._array[slice_],
                              unit=self._unit,
                              parent=self._parent,
                              name=self._name)

    def __len__(self):
        if self.shape:
            return self.shape[0]
        else:
            return 0

    def __str__(self):
        return "Vector<{}, (x,y,z)>".format(self._make_string())

    @property
    def values(self):
        return self._array

    @property
    def norm(self):
        return Array(values=np.linalg.norm(self._array, axis=-1), unit=self.unit)

    # @property
    # def nvec(self):
    #     if self._array.shape:
    #         if len(self._array.shape) == 2:
    #             return self._array.shape[-1]
    #         else:
    #             return 1
    #     return 0

    @property
    def shape(self):
        return self._array[..., 0].shape

    @property
    def x(self):
        return Array(values=self._array[..., 0],
                     unit=self._unit,
                     parent=self._parent,
                     name=self._name + "_x")

    @property
    def y(self):
        if self._array.shape[-1] > 1:
            return Array(values=self._array[..., 1],
                         unit=self._unit,
                         parent=self._parent,
                         name=self._name + "_y")

    @property
    def z(self):
        if self._array.shape[-1] > 2:
            return Array(values=self._array[..., 2],
                         unit=self._unit,
                         parent=self._parent,
                         name=self._name + "_z")

    def _binary_op(self, op, lhs, rhs, mul_or_div=False, out=None):
        if isinstance(rhs, Quantity):
            rhs = Array(values=rhs.magnitude, unit=1.0 * rhs.units)
        if isinstance(rhs, (int, float, np.ndarray)):
            rhs = Array(values=rhs)
        if not hasattr(rhs, "x"):
            rhs = rhs.reshape(rhs.shape + (1, ))
        return self._compute_binary_op(op, lhs, rhs, mul_or_div, out)

    def __radd__(self, other):
        return add(self, other)

    def __rsub__(self, other):
        return -sub(self, other)

    def min(self):
        return Array(values=self.norm._array.min(), unit=self._unit)

    def max(self):
        return Array(values=self.norm._array.max(), unit=self._unit)
