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


class Vector:
    def __init__(self, x, y=None, z=None, parent=None, name=""):

        self._parent = parent
        self._name = name

        self.x = x
        self.x.name = self._name + "_x"

        self.y = None
        if y is not None:
            self.y = y
            self.y.name = self._name + "_y"
            if self.y.shape != self.x.shape:
                raise ValueError("The shape of component y must be the same as the "
                                 "shape of the x component")

        self.z = None
        if z is not None:
            self.z = z
            self.z.name = self._name + "_z"
            if self.z.shape != self.x.shape:
                raise ValueError("The shape of component z must be the same as the "
                                 "shape of the x component")

    def __getitem__(self, slice_):
        slice_ = tuple((slice_, )) + (slice(None, None, None), )
        x = self.x[slice_]
        y = self.y[slice_] if self.y is not None else None
        z = self.z[slice_] if self.z is not None else None
        return self.__class__(x=x,
                              y=y,
                              z=z,
                              unit=self._unit,
                              parent=self._parent,
                              name=self._name)

    def __len__(self):
        return len(self.x)

    def __str__(self):
        name_str = "'" + self._name + "' "
        if len(self) == 0:
            values_str = "Value: " + value_to_string(self.values)
        else:
            values_str = "Min: " + value_to_string(
                self.norm.min().values) + " Max: " + value_to_string(
                    self.norm.max().values)
        unit_str = " [{:~}] ".format(self.unit.units)
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
        x = self.x.copy()
        y = self.y.copy() if self.y is not None else None
        z = self.z.copy() if self.z is not None else None
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
        if (self.y is None) and (self.z is None):
            return self.x
        comps = [self.x.array, self.y.array]
        if self.z is not None:
            comps.append(self.z.array)
        return Array(values=np.linalg.norm(np.asarray(comps).T, axis=-1),
                     unit=self.unit)

    @property
    def unit(self):
        return self.x.unit

    @unit.setter
    def unit(self, unit_):
        self.x.unit = unit_
        if self.y is not None:
            self.y = unit_
        if self.z is not None:
            self.z = unit_

    @property
    def ndim(self):
        return self.x.ndim

    @property
    def nvec(self):
        if (self.y is None) and (self.z is None):
            return 1
        if self.z is None:
            return 2
        return 3

    @property
    def shape(self):
        return self.x.shape

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent_):
        self._parent = parent_

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_):
        self._name = name_

    @property
    def label(self):
        return make_label(name=self._name, unit=self._unit.units)

    def _to_vector(self, rhs):
        if isinstance(rhs, Quantity):
            rhs = Array(values=rhs.magnitude, unit=1.0 * rhs.units)
        if isinstance(rhs, (int, float, np.ndarray)):
            rhs = Array(values=rhs)
        if isinstance(rhs, Array):
            rhs = self.__class__(**{
                c: rhs
                for c in "xyz" if getattr(self, c) is not None
            },
                                 unit=rhs.unit)
        if self.nvec != rhs.nvec:
            raise ValueError("Operands do not have the same number of components.")
        return rhs
        # if isinstance(rhs, Quantity):
        #     rhs = Array(values=rhs.magnitude, unit=1.0 * rhs.units)
        # if isinstance(rhs, (int, float, np.ndarray)):
        #     rhs = Array(values=rhs)
        # if not hasattr(rhs, "x"):
        #     rhs = rhs.reshape(rhs.shape + (1, ))
        # return self._compute_binary_op(op, lhs, rhs, mul_or_div, out)

    def __add__(self, other):
        other = self._to_vector(other)
        x = (self.x + other.x)
        y = (self.y + other.y) if self.y is not None else None
        z = (self.z + other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __iadd__(self, other):
        other = self._to_vector(other)
        self.x += other.x
        if self.y is not None:
            self.y += other.y
        if self.z is not None:
            self.z += other.z

    def __sub__(self, other):
        other = self._to_vector(other)
        x = (self.x - other.x)
        y = (self.y - other.y) if self.y is not None else None
        z = (self.z - other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __isub__(self, other):
        other = self._to_vector(other)
        self.x -= other.x
        if self.y is not None:
            self.y -= other.y
        if self.z is not None:
            self.z -= other.z

    def __mul__(self, other):
        other = self._to_vector(other)
        x = (self.x * other.x)
        y = (self.y * other.y) if self.y is not None else None
        z = (self.z * other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __imul__(self, other):
        other = self._to_vector(other)
        self.x *= other.x
        if self.y is not None:
            self.y *= other.y
        if self.z is not None:
            self.z *= other.z

    def __truediv__(self, other):
        other = self._to_vector(other)
        x = (self.x / other.x)
        y = (self.y / other.y) if self.y is not None else None
        z = (self.z / other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __itruediv__(self, other):
        other = self._to_vector(other)
        self.x /= other.x
        if self.y is not None:
            self.y /= other.y
        if self.z is not None:
            self.z /= other.z

    # def __rmul__(self, other):
    #     return self * other

    # def __rtruediv__(self, other):
    #     out = np.reciprocal(self / other)
    #     out._unit = 1.0 / out._unit
    #     return out

    # def __radd__(self, other):
    #     return add(self, other)

    # def __rsub__(self, other):
    #     return -sub(self, other)

    def __pow__(self, number):
        x = self.x**number
        y = self.y**number if self.y is not None else None
        z = self.z**number if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __neg__(self):
        x = -self.x
        y = -self.y if self.y is not None else None
        z = -self.z if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __lt__(self, other):
        other = self._to_vector(other)
        x = (self.x < other.x)
        y = (self.y < other.y) if self.y is not None else None
        z = (self.z < other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __le__(self, other):
        other = self._to_vector(other)
        x = (self.x <= other.x)
        y = (self.y <= other.y) if self.y is not None else None
        z = (self.z <= other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __gt__(self, other):
        other = self._to_vector(other)
        x = (self.x > other.x)
        y = (self.y > other.y) if self.y is not None else None
        z = (self.z > other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __ge__(self, other):
        other = self._to_vector(other)
        x = (self.x >= other.x)
        y = (self.y >= other.y) if self.y is not None else None
        z = (self.z >= other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __eq__(self, other):
        other = self._to_vector(other)
        x = (self.x == other.x)
        y = (self.y == other.y) if self.y is not None else None
        z = (self.z == other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def __ne__(self, other):
        other = self._to_vector(other)
        x = (self.x != other.x)
        y = (self.y != other.y) if self.y is not None else None
        z = (self.z != other.z) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def to(self, unit):
        x = self.x.to(unit)
        y = self.y.to(unit) if self.y is not None else None
        z = self.z.to(unit) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

    def _wrap_numpy(self, func, *args, **kwargs):
        # # if func.__name__ in self.special_functions:
        # #     unit = func(self.unit, *args[1:], **kwargs)
        # # else:
        # #     unit = self.unit
        # if isinstance(args[0], tuple) or isinstance(args[0], list):
        #     # Case where we have a sequence of arrays, e.g. `concatenate`
        #     for a in args[0]:
        #         if a.unit != unit:
        #             raise TypeError("Could not {} types {} and {}.".format(
        #                 func.__name__, self, a))
        #     args = (tuple(a._array for a in args[0]), ) + args[1:]
        # elif (len(args) > 1 and hasattr(args[1], "_array")):
        #     if hasattr(args[0], "_array"):
        #         # Case of a binary operation, with two Arrays, e.g. `dot`
        #         # TODO: what should we do with the unit? Apply the func to it?
        #         unit = func(args[0].unit, args[1].unit, *args[2:], **kwargs)
        #         args = (args[0]._array, args[1]._array) + args[2:]
        #     else:
        #         # Case of a binary operation: ndarray with Array
        #         # In this case, only multiply is allowed?
        #         if func.__name__ != "multiply":
        #             raise RuntimeError("Cannot use operation {} between ndarray and "
        #                                "Array".format(func.__name__))
        #         args = (args[0], args[1]._array) + args[2:]
        # else:
        #     args = (args[0]._array, ) + args[1:]
        x = func(self.x, *args[1:], **kwargs)
        y = func(self.y, *args[1:], **kwargs) if self.y is not None else None
        z = func(self.z, *args[1:], **kwargs) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)

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

    def min(self):
        return Array(values=self.norm._array.min(), unit=self._unit)

    def max(self):
        return Array(values=self.norm._array.max(), unit=self._unit)

    def reshape(self, *shape):
        x = self.x.reshape(*shape)
        y = self.y.reshape(*shape) if self.y is not None else None
        z = self.z.reshape(*shape) if self.z is not None else None
        return self.__class__(x=x, y=y, z=z)
