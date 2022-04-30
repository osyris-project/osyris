# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
# from .operators import add, iadd, sub, isub, mul, imul, div, idiv, comp
from .tools import value_to_string, make_label
from .. import units


class Base:
    def __init__(self, unit=None, parent=None, name=""):

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

    def __repr__(self):
        return str(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, unit_):
        self._unit = unit_

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

    def _compute_binary_op(self, op, lhs, rhs, mul_or_div, out):
        ratio = op(lhs._unit, rhs._unit) if mul_or_div else rhs.unit.to(lhs._unit.units)
        result = op(lhs._array, rhs._array, out=out)
        result *= ratio.magnitude
        if out is None:
            return lhs.__class__(values=result, unit=1.0 * ratio.units)
        else:
            if mul_or_div:
                lhs._unit = 1.0 * ratio.units
            return lhs

    def _binary_op(self, op, lhs, rhs, mul_or_div=False, out=None):
        if isinstance(rhs, Quantity):
            rhs = lhs.__class__(values=rhs.magnitude, unit=1.0 * rhs.units)
        if isinstance(rhs, (int, float, np.ndarray)):
            rhs = lhs.__class__(values=rhs)
        if (len(rhs) > 1) and (lhs.ndim != rhs.ndim):
            return NotImplemented
        return self._compute_binary_op(op, lhs, rhs, mul_or_div, out)

    def __add__(self, other):
        return self._binary_op(np.add, self, other)

    def __iadd__(self, other):
        return self._binary_op(np.add, self, other, out=self._array)

    def __sub__(self, other):
        return self._binary_op(np.subtract, self, other)

    def __isub__(self, other):
        return self._binary_op(np.subtract, self, other, out=self._array)

    def __mul__(self, other):
        return self._binary_op(np.multiply, self, other, mul_or_div=True)

    def __imul__(self, other):
        return self._binary_op(np.multiply,
                               self,
                               other,
                               mul_or_div=True,
                               out=self._array)

    def __truediv__(self, other):
        return self._binary_op(np.divide, self, other, mul_or_div=True)

    def __itruediv__(self, other):
        return self._binary_op(np.divide, self, other, mul_or_div=True, out=self._array)

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        out = np.reciprocal(self / other)
        out._unit = 1.0 / out._unit
        return out

    def __pow__(self, number):
        return np.power(self, number)

    def __neg__(self):
        return np.negative(self)

    def __lt__(self, other):
        return _comparison_op(np.less, self, other)

    def __le__(self, other):
        return _comparison_op(np.less_equal, self, other)

    def __gt__(self, other):
        return _comparison_op(np.greater, self, other)

    def __ge__(self, other):
        return _comparison_op(np.greater_equal, self, other)

    def __eq__(self, other):
        return _comparison_op(np.equal, self, other)

    def __ne__(self, other):
        return _comparison_op(np.not_equal, self, other)

    def to(self, unit):
        if isinstance(unit, str):
            new_unit = units(unit)
        else:
            new_unit = unit
        ratio = self._unit.to(new_unit) / new_unit
        return self.__class__(values=self._array * ratio.magnitude, unit=1.0 * new_unit)

    def _wrap_numpy(self, func, *args, **kwargs):
        if func.__name__ in self.special_functions:
            unit = func(self.unit, *args[1:], **kwargs)
        else:
            unit = self.unit
        if isinstance(args[0], tuple) or isinstance(args[0], list):
            # Case where we have a sequence of arrays, e.g. `concatenate`
            for a in args[0]:
                if a.unit != unit:
                    raise TypeError("Could not {} types {} and {}.".format(
                        func.__name__, self, a))
            args = (tuple(a._array for a in args[0]), ) + args[1:]
        elif (len(args) > 1 and hasattr(args[1], "_array")):
            if hasattr(args[0], "_array"):
                # Case of a binary operation, with two Arrays, e.g. `dot`
                # TODO: what should we do with the unit? Apply the func to it?
                unit = func(args[0].unit, args[1].unit, *args[2:], **kwargs)
                args = (args[0]._array, args[1]._array) + args[2:]
            else:
                # Case of a binary operation: ndarray with Array
                # In this case, only multiply is allowed?
                if func.__name__ != "multiply":
                    raise RuntimeError("Cannot use operation {} between ndarray and "
                                       "Array".format(func.__name__))
                args = (args[0], args[1]._array) + args[2:]
        else:
            args = (args[0]._array, ) + args[1:]
        result = func(*args, **kwargs)
        return self.__class__(values=result, unit=unit)

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
        return self.__class__(values=self._array.min(), unit=self._unit)

    def max(self):
        return self.__class__(values=self._array.max(), unit=self._unit)

    def reshape(self, *shape):
        return self.__class__(values=self._array.reshape(*shape), unit=self._unit)
