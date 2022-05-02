# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
from .tools import value_to_string, make_label
from .. import units


def binary_op(op, lhs, rhs, mul_or_div=False, out=None):
    if isinstance(rhs, Quantity):
        rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
    if isinstance(rhs, (int, float, np.ndarray)):
        rhs = lhs.__class__(values=rhs)
    if not isinstance(rhs, lhs.__class__):
        return NotImplemented

    ratio = op(lhs.unit, rhs.unit) if mul_or_div else rhs.unit.to(lhs.unit.units)
    result = op(lhs._array, rhs._array, out=out)
    result *= ratio.magnitude
    if out is None:
        return lhs.__class__(values=result, unit=1.0 * ratio.units)
    else:
        if mul_or_div:
            lhs.unit = 1.0 * ratio.units
        return lhs


class Array:
    def __init__(self, values=0, unit=None, parent=None, name=""):

        self._array = np.asarray(values)

        if unit is None:
            self.unit = 1.0 * units.dimensionless
        elif isinstance(unit, str):
            self.unit = units(unit)
        elif isinstance(unit, Quantity):
            self.unit = unit
        elif isinstance(unit, Unit):
            self.unit = 1.0 * unit
        else:
            raise TypeError("Unsupported unit type {}".format(type(unit)))
        self.parent = parent
        self.name = name
        self._special_functions = ["sqrt", "power", "multiply", "divide"]

    def __getitem__(self, slice_):
        return self.__class__(values=self._array[slice_],
                              unit=self.unit,
                              parent=self.parent,
                              name=self.name)

    def __len__(self):
        if self._array.shape:
            return len(self._array)
        else:
            return 0

    def __str__(self):
        name_str = "'" + self.name + "' "
        if len(self) == 0:
            values_str = "Value: " + value_to_string(self.values)
        else:
            values_str = "Min: " + value_to_string(
                self.min().values) + " Max: " + value_to_string(self.max().values)
        unit_str = " [{:~}] ".format(self.unit.units)
        shape_str = str(self.shape)
        return "Array<" + name_str + values_str + unit_str + shape_str + ">"

    def __repr__(self):
        return str(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def copy(self):
        return self.__class__(values=self._array.copy(),
                              unit=self.unit.copy(),
                              name=str(self.name))

    @property
    def values(self):
        if not self._array.shape:
            return self._array[()]
        else:
            return self._array

    @values.setter
    def values(self, values_):
        self._array = values_

    @property
    def norm(self):
        if self._array.ndim < 2:
            return self
        else:
            return self.__class__(values=np.linalg.norm(self._array, axis=1),
                                  unit=self.unit)

    @property
    def ndim(self):
        return self._array.ndim

    @property
    def shape(self):
        return self._array.shape

    @property
    def label(self):
        return make_label(name=self.name, unit=self.unit.units)

    def _to_array(self, other):
        if isinstance(other, Quantity):
            other = self.__class__(values=other.magnitude, unit=1.0 * other.units)
        if isinstance(other, self.__class__):
            if other.unit != self.unit:
                other = other.to(self.unit)
            return other

    def __add__(self, other):
        return binary_op(np.add, self, other)

    def __iadd__(self, other):
        return binary_op(np.add, self, other, out=self._array)

    def __sub__(self, other):
        return binary_op(np.subtract, self, other)

    def __isub__(self, other):
        return binary_op(np.subtract, self, other, out=self._array)

    def __mul__(self, other):
        return binary_op(np.multiply, self, other, mul_or_div=True)

    def __imul__(self, other):
        return binary_op(np.multiply, self, other, mul_or_div=True, out=self._array)

    def __truediv__(self, other):
        return binary_op(np.divide, self, other, mul_or_div=True)

    def __itruediv__(self, other):
        return binary_op(np.divide, self, other, mul_or_div=True, out=self._array)

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        out = np.reciprocal(self / other)
        out.unit = 1.0 / out.unit
        return out

    def __pow__(self, number):
        return np.power(self, number)

    def __neg__(self):
        return np.negative(self)

    def __lt__(self, other):
        return np.less(self, self._to_array(other))

    def __le__(self, other):
        return np.less_equal(self, self._to_array(other))

    def __gt__(self, other):
        return np.greater(self, self._to_array(other))

    def __ge__(self, other):
        return np.greater_equal(self, self._to_array(other))

    def __eq__(self, other):
        return np.equal(self, self._to_array(other))

    def __ne__(self, other):
        return np.not_equal(self, self._to_array(other))

    def to(self, unit):
        if isinstance(unit, str):
            new_unit = units(unit)
        else:
            new_unit = unit
        ratio = self.unit.to(new_unit) / new_unit
        return self.__class__(values=self._array * ratio.magnitude, unit=1.0 * new_unit)

    def _extract_underlying(self, args):
        return tuple(a._array if isinstance(a, self.__class__) else a for a in args)

    def _wrap_numpy(self, func, *args, **kwargs):
        if func.__name__ in self._special_functions:
            unit = func(self.unit, *args[1:], **kwargs)
        else:
            unit = self.unit
        if isinstance(args[0], (tuple, list)):
            args = (self._extract_underlying(args[0]), ) + self._extract_underlying(
                args[1:])
        else:
            args = self._extract_underlying(args)
        result = func(*args, **kwargs)
        return self.__class__(values=result,
                              unit=unit if result.dtype in (int, float) else None)

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
        return self.__class__(values=self._array.min(), unit=self.unit)

    def max(self):
        return self.__class__(values=self._array.max(), unit=self.unit)

    def reshape(self, *shape):
        return self.__class__(values=self._array.reshape(*shape), unit=self.unit)

    def nbytes(self):
        return self._array.nbytes