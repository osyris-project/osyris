# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
from .array import Array
from .operators import add, sub
from .tools import value_to_string, make_label
from .. import units


class Vector:
    def __init__(self, x, y=None, z=None, unit=None, parent=None, name=""):

        self._x = np.asarray(x)
        self._y = np.asarray(y) if y is not None else None
        self._z = np.asarray(z) if z is not None else None

        if unit is None:
            self._unit = 1.0 * units.dimensionless
        elif isinstance(unit, str):
            self._unit = units(unit)
        elif isinstance(unit, Quantity):
            self._unit = unit
        elif isinstance(unit, Unit):
            self._unit = 1.0 * unit
        else:
            raise TypeError("Unsupported unit type {}".format(type(unit)))
        self._parent = parent
        self._name = name
        # self.special_functions = ["sqrt", "power"]

    def __getitem__(self, slice_):
        slice_ = tuple((slice_, )) + (slice(None, None, None), )
        x = self._x[slice_]
        y = self._y[slice_] if self._y is not None else None
        z = self._z[slice_] if self._z is not None else None

        return self.__class__(x=x,
                              y=y,
                              z=z,
                              unit=self._unit,
                              parent=self._parent,
                              name=self._name)

    def __len__(self):
        if self._x.shape:
            return self._x.shape[0]
        else:
            return 0

    def __str__(self):
        name_str = "'" + self._name + "' "
        if len(self) == 0:
            values_str = "Value: " + value_to_string(self.values)
        else:
            values_str = "Min: " + value_to_string(
                self.min().values) + " Max: " + value_to_string(self.max().values)
        unit_str = " [{:~}] ".format(self._unit.units)
        shape_str = str(self.shape)
        return "Vector<" + name_str + values_str + unit_str + shape_str

        # return "Vector<{}, (x,y,z)>".format(self._make_string())

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
