# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint import Quantity
from pint.errors import DimensionalityError
from .base import Base
from .tools import value_to_string
from .. import units

APPLY_OP_TO_UNIT = ("multiply", "true_divide", "divide", "sqrt", "power", "reciprocal")


def _binary_op(op, lhs, rhs, strict=True, **kwargs):
    if not isinstance(rhs, lhs.__class__):
        try:
            rhs = lhs.__class__(rhs)
        except NotImplementedError:
            return NotImplemented
    if strict:
        rhs = rhs.to(lhs.unit)
    else:
        try:
            rhs = rhs.to(lhs.unit)
        except DimensionalityError:
            pass
    return op(lhs, rhs, **kwargs)


class Array(Base):
    def __init__(self, values, unit=None, parent=None, name=""):
        if isinstance(values, Base):
            raise NotImplementedError("Cannot create Array from Array or Vector.")

        if isinstance(values, Quantity):
            if unit is not None:
                raise ValueError(
                    "Cannot set unit when creating an Array from a Quantity."
                )
            self._array = values.magnitude
            self.unit = values.units
        else:
            self._array = values
            self.unit = units(unit)
        if not isinstance(self._array, np.ndarray):
            self._array = np.asarray(self._array)

        self.parent = parent
        self.name = name

    def __getitem__(self, slice_):
        if isinstance(slice_, Base):
            if not isinstance(slice_, self.__class__):
                raise ValueError(
                    "Cannot slice using a Vector, only Array is supported."
                )
            if slice_.dtype not in ("int32", "int64", bool):
                raise TypeError(
                    "Dtype of the Array must be integer or bool for slicing."
                )
            slice_ = slice_.values
        return self.__class__(
            values=self._array[slice_],
            unit=self.unit,
            parent=self.parent,
            name=self.name,
        )

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
            values_str = (
                "Min: "
                + value_to_string(self.min().values)
                + " Max: "
                + value_to_string(self.max().values)
            )
        unit_str = " [{:~}] ".format(self.unit)
        shape_str = str(self.shape)
        return name_str + values_str + unit_str + shape_str

    def copy(self):
        return self.__class__(
            values=self._array.copy(), unit=units(self.unit), name=str(self.name)
        )

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
        return self

    @property
    def ndim(self):
        return self._array.ndim

    @property
    def shape(self):
        return self._array.shape

    @property
    def dtype(self):
        return self._array.dtype

    def __add__(self, other):
        return _binary_op(np.add, self, other)

    def __iadd__(self, other):
        return _binary_op(np.add, self, other, out=self)

    def __sub__(self, other):
        return _binary_op(np.subtract, self, other)

    def __isub__(self, other):
        return _binary_op(np.subtract, self, other, out=self)

    def __mul__(self, other):
        return _binary_op(np.multiply, self, other, strict=False)

    def __imul__(self, other):
        return _binary_op(np.multiply, self, other, strict=False, out=self)

    def __truediv__(self, other):
        return _binary_op(np.divide, self, other, strict=False)

    def __itruediv__(self, other):
        return _binary_op(np.divide, self, other, strict=False, out=self)

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        return np.reciprocal(self / other)

    def __pow__(self, number):
        return np.power(self, number)

    def __neg__(self):
        return np.negative(self)

    def __lt__(self, other):
        return _binary_op(np.less, self, other)

    def __le__(self, other):
        return _binary_op(np.less_equal, self, other)

    def __gt__(self, other):
        return _binary_op(np.greater, self, other)

    def __ge__(self, other):
        return _binary_op(np.greater_equal, self, other)

    def __eq__(self, other):
        return _binary_op(np.equal, self, other)

    def __ne__(self, other):
        return _binary_op(np.not_equal, self, other)

    def __and__(self, other):
        return _binary_op(np.logical_and, self, other)

    def __or__(self, other):
        return _binary_op(np.logical_or, self, other)

    def __xor__(self, other):
        return _binary_op(np.logical_xor, self, other)

    def __invert__(self):
        return np.logical_not(self)

    def to(self, unit):
        new_unit = units(unit)
        if self.unit == new_unit:
            return self
        ratio = (1.0 * self.unit).to(new_unit) / (1.0 * new_unit)
        return self.__class__(values=self._array * ratio.magnitude, unit=new_unit)

    def _maybe_array(self, arg):
        if isinstance(arg, self.__class__):
            return arg._array
        if isinstance(arg, Quantity):
            return arg.magnitude
        return arg

    def _extract_arrays_from_args(self, args):
        return tuple(self._maybe_array(a) for a in args)

    def _extract_arrays_from_kwargs(self, kwargs):
        return {key: self._extract_arrays_from_args(a) for key, a in kwargs.items()}

    def _maybe_unit(self, arg):
        if hasattr(arg, "unit"):
            return 1.0 * arg.unit
        if hasattr(arg, "units"):
            return 1.0 * arg.units
        return arg

    def _extract_units(self, args):
        return tuple(self._maybe_unit(a) for a in args)

    def _wrap_numpy(self, func, *args, **kwargs):
        if isinstance(args[0], (tuple, list)):
            array_args = (
                self._extract_arrays_from_args(args[0]),
            ) + self._extract_arrays_from_args(args[1:])
        else:
            array_args = self._extract_arrays_from_args(args)
        result = func(*array_args, **self._extract_arrays_from_kwargs(kwargs))

        unit = None
        if result.dtype in (int, float):
            if func.__name__ in APPLY_OP_TO_UNIT:
                unit = func(
                    *self._extract_units(args),
                    **{key: a for key, a in kwargs.items() if key != "out"},
                ).units
            else:
                unit = self.unit

        if "out" in kwargs:
            kwargs["out"][0].unit = unit
            return kwargs["out"][0]
        else:
            return self.__class__(values=result, unit=unit)

    def reshape(self, *shape):
        return self.__class__(values=self._array.reshape(*shape), unit=self.unit)

    @property
    def nbytes(self):
        return self._array.nbytes
