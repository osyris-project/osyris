# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
from .array import Array
# from .base import Base
# from .operators import add, sub
from .tools import value_to_string, make_label
from .. import units


class Vector:
    def __init__(self, x, y=None, z=None, parent=None, name=""):
        self.parent = parent
        unit = x.unit
        self.x = Array(values=x.values, unit=unit)
        self.y = self._validate_component(y, unit=unit)
        self.z = self._validate_component(z, unit=unit)
        self.name = name

    @classmethod
    def from_values(cls, values, unit=None):
        assert values.ndim > 1
        nvec = values.shape[-1]
        x = Array(values=values[..., 0], unit=unit)
        y = Array(values=values[..., 1], unit=unit) if nvec > 1 else None
        z = Array(values=values[..., 2], unit=unit) if nvec > 2 else None
        return cls(x=x, y=y, z=z)

    def _validate_component(self, array, unit):
        if array is None:
            return array
        if array.shape != self.x.shape:
            raise ValueError(f"The shape of the component does not match the "
                             "shape of the x component")
        if array.unit != unit:
            raise ValueError(f"The unit of the component does not match the "
                             "unit of the x component")
        return Array(values=array.values, unit=unit)

    @property
    def _xyz(self):
        out = {'x': self.x}
        if self.y is not None:
            out['y'] = self.y
        if self.z is not None:
            out['z'] = self.z
        return out

    def __getitem__(self, slice_):
        return self.__class__(**{c: xyz[slice_]
                                 for c, xyz in self._xyz.items()},
                              parent=self.parent,
                              name=self._name)

    def __len__(self):
        return len(self.x)

    def __str__(self):
        name_str = "'" + self._name + "' "
        xyz = self._xyz
        if len(self) == 0:
            values_str = "Value: " + ",".join(
                value_to_string(x.values) for x in xyz.values())
        else:
            norm = self.norm
            values_str = "Min: " + value_to_string(
                norm.min().values) + " Max: " + value_to_string(norm.max().values)
        unit_str = " [{:~}] ".format(self.unit.units)
        shape_str = str(self.shape)
        comps_str = ", (" + ",".join(x for x in xyz) + ")>"
        return "Vector<" + name_str + values_str + unit_str + shape_str + comps_str

    def __repr__(self):
        return str(self)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def copy(self):
        return self.__class__(**{c: xyz.copy()
                                 for c, xyz in self._xyz.items()},
                              name=str(self._name))

    @property
    def norm(self):
        if (self.y is None) and (self.z is None):
            return self.x
        out = self.x.values * self.x.values
        out += self.y.values * self.y.values
        if self.z is not None:
            out += self.z.values * self.z.values
        return Array(values=np.sqrt(out), unit=self.x.unit)

    @property
    def unit(self):
        return self.x.unit

    @unit.setter
    def unit(self, unit_):
        self.x.unit = unit_
        if self.y is not None:
            self.y.unit = unit_
        if self.z is not None:
            self.z.unit = unit_

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
    def name(self):
        return self._name

    @name.setter
    def name(self, name_):
        self._name = name_
        for c, xyz in self._xyz.items():
            xyz.name = self._name + "_" + c

    @property
    def label(self):
        return make_label(name=self._name, unit=self.unit.units)

    def _to_vector(self, rhs):
        if isinstance(rhs, Quantity):
            rhs = Array(values=rhs.magnitude, unit=1.0 * rhs.units)
        if isinstance(rhs, (int, float, np.ndarray)):
            rhs = Array(values=rhs)
        if isinstance(rhs, Array):
            rhs = self.__class__(**{c: rhs for c in self._xyz.keys()})
        if self.nvec != rhs.nvec:
            raise ValueError("Operands do not have the same number of components.")
        return rhs

    def __add__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz + getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __iadd__(self, other):
        other = self._to_vector(other)
        self.x += other.x
        if self.y is not None:
            self.y += other.y
        if self.z is not None:
            self.z += other.z
        return self

    def __sub__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz - getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __isub__(self, other):
        other = self._to_vector(other)
        self.x -= other.x
        if self.y is not None:
            self.y -= other.y
        if self.z is not None:
            self.z -= other.z
        return self

    def __mul__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz * getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __imul__(self, other):
        other = self._to_vector(other)
        self.x *= other.x
        if self.y is not None:
            self.y *= other.y
        if self.z is not None:
            self.z *= other.z
        return self

    def __truediv__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz / getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __itruediv__(self, other):
        other = self._to_vector(other)
        self.x /= other.x
        if self.y is not None:
            self.y /= other.y
        if self.z is not None:
            self.z /= other.z
        return self

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        out = np.reciprocal(self / other)
        out.unit = 1.0 / out.unit
        return out

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return -(self - other)

    def __pow__(self, number):
        return self.__class__(**{c: xyz**number for c, xyz in self._xyz.items()})

    def __neg__(self):
        return self.__class__(**{c: -xyz for c, xyz in self._xyz.items()})

    def __lt__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz < getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __le__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz <= getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __gt__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz > getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __ge__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz >= getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __eq__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz == getattr(other, c)
               for c, xyz in self._xyz.items()})

    def __ne__(self, other):
        other = self._to_vector(other)
        return self.__class__(
            **{c: xyz != getattr(other, c)
               for c, xyz in self._xyz.items()})

    def to(self, unit):
        return self.__class__(**{c: xyz.to(unit) for c, xyz in self._xyz.items()})

    def _wrap_numpy(self, func, *args, **kwargs):
        if isinstance(args[0], (tuple, list)):
            # # Case where we have a sequence of vectors, e.g. `concatenate`
            out = {
                c: func(tuple(getattr(a, c) for a in args[0]), *args[1:], **kwargs)
                for c, xyz in args[0][0]._xyz.items()
            }
        elif (len(args) > 1 and isinstance(args[1], self.__class__)):
            # Case of a binary operation, with two vectors, e.g. `dot`
            out = {
                c: func(xyz, getattr(args[1], c), *args[2:], **kwargs)
                for c, xyz in args[0]._xyz.items()
            }
        else:
            out = {c: func(xyz, *args[1:], **kwargs) for c, xyz in args[0]._xyz.items()}
        return self.__class__(**out)

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
        return np.amin(self)

    def max(self):
        return np.amax(self)

    def reshape(self, *shape):
        return self.__class__(
            **{c: xyz.reshape(*shape)
               for c, xyz in self._xyz.items()})

    @property
    def nbytes(self):
        return np.sum([xyz.nbytes for xyz in self._xyz.values()])

    def dot(self, other):
        out = np.zeros(self.shape)
        for (c1, c2) in zip(self._xyz.values(), other._xyz.values()):
            out += (c1 * c2).values
        return Array(values=out, unit=self.unit * other.unit)
