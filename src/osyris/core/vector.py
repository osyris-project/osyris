# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
from .array import Array
from .base import Base
from .operators import add, sub
from .tools import value_to_string, make_label
from .. import units


class Vector(Base):
    def __init__(self, x, y=None, z=None, unit=None, parent=None, name=""):
        super().__init__(unit=unit, parent=parent, name=name)

        self._x = x
        self._x.name = self._name + "_x"

        self._y = None
        if y is not None:
            self._y = y
            self._y.name = self._name + "_y"
            if self._y.shape != self._x.shape:
                raise ValueError("The shape of component y must be the same as the "
                                 "shape of the x component")

        self._z = None
        if z is not None:
            self._z = z
            self._z.name = self._name + "_z"
            if self._z.shape != self._x.shape:
                raise ValueError("The shape of component z must be the same as the "
                                 "shape of the x component")

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
        return len(self._x)

    def __str__(self):
        name_str = "'" + self._name + "' "
        if len(self) == 0:
            values_str = "Value: " + value_to_string(self.values)
        else:
            values_str = "Min: " + value_to_string(
                self.norm.min().values) + " Max: " + value_to_string(
                    self.norm.max().values)
        unit_str = " [{:~}] ".format(self._unit.units)
        shape_str = str(self.shape)
        return "Vector<" + name_str + values_str + unit_str + shape_str + ", (x,y,z)>"

        # return "Vector<{}, (x,y,z)>".format(self._make_string())

    def __repr__(self):
        return str(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def copy(self):
        x = self._x.copy()
        y = self._y.copy() if self._y is not None else None
        z = self._z.copy() if self._z is not None else None
        return self.__class__(x=x,
                              y=y,
                              z=z,
                              unit=self._unit.copy(),
                              name=str(self._name))

    # @property
    # def values(self):
    #     return self._array

    @property
    def norm(self):
        if (self._y is None) and (self._z is None):
            return self._x
        comps = [self._x, self._y]
        if self._z is not None:
            comps.append(self._z)
        return Array(values=np.linalg.norm(np.asarray(comps).T, axis=-1),
                     unit=self.unit)


    @property
    def ndim(self):
        return self._x.ndim

    @property
    def nvec(self):
        if (self._y is None) and (self._z is None):
            return 1
        if self._z is None:
            return 2
        return 3

    @property
    def shape(self):
        return self._x.shape


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_):
        self._name = name_

    @property
    def label(self):
        return make_label(name=self._name, unit=self._unit.units)

    @property
    def x(self):
        return Array(values=self._x,
                     unit=self._unit,
                     parent=self._parent,
                     name=self._name + "_x")

    @property
    def y(self):
        if self._y is not None:
            return Array(values=self._y,
                         unit=self._unit,
                         parent=self._parent,
                         name=self._name + "_y")

    @property
    def z(self):
        if self._z is not None:
            return Array(values=self._z,
                         unit=self._unit,
                         parent=self._parent,
                         name=self._name + "_z")

    def _binary_op(self, op, lhs, rhs, mul_or_div=False, out=None):
        x = op(


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
