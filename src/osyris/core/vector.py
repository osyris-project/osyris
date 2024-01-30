# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
import numpy as np
from pint import Quantity
from .base import Base
from .array import Array
from .tools import value_to_string


def _binary_op(op, lhs, rhs):
    if isinstance(rhs, (int, float, np.ndarray, Quantity)):
        rhs = Array(values=rhs)
    if isinstance(rhs, Array):
        rhs = lhs.__class__(**{c: rhs for c in lhs._xyz.keys()})
    if lhs.nvec != rhs.nvec:
        raise ValueError("Operands do not have the same number of components.")

    return lhs.__class__(
        **{c: getattr(xyz, op)(getattr(rhs, c))
           for c, xyz in lhs._xyz.items()})


def _get_colatitude(pos):
    colatitude = np.arctan2(pos.cyl_r, pos.z).values
    return Array(values=colatitude, unit="rad", name=pos.name + "_theta")


def _get_azimuth(pos):
    azimuth = np.arctan2(pos.y.values, pos.x.values)
    return Array(values=azimuth, unit="rad", name=pos.name + "_phi")


class Vector(Base):

    def __init__(self, x, y=None, z=None, parent=None, name="", unit=None):

        if isinstance(x, Array):
            if unit is not None:
                raise ValueError("Can only set unit when creating Vector from values.")
            unit = x.unit
            x = x.values
            y = self._validate_component(y, x.shape, unit)
            z = self._validate_component(z, x.shape, unit)
        self.x = Array(values=x, unit=unit)
        self.y = Array(values=y, unit=unit) if y is not None else None
        self.z = Array(values=z, unit=unit) if z is not None else None

        self.parent = parent
        self.name = name

    def _validate_component(self, array, shape, unit):
        if array is None:
            return
        if array.shape != shape:
            raise ValueError("The shape of the component does not match the "
                             "shape of the x component")
        if array.unit != unit:
            raise ValueError("The unit of the component does not match the "
                             "unit of the x component")
        return array.values

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
        comps_str = ", {" + ",".join(x for x in xyz) + "}"
        if len(self) == 0:
            values_str = "Value: " + ", ".join(
                value_to_string(x.values) for x in xyz.values())
            unit_str = " [{:~}] ".format(self.unit)
            shape_str = str(self.shape)
            return name_str + values_str + unit_str + shape_str + comps_str
        else:
            return str(self.norm) + comps_str

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
        return Array(values=np.sqrt(out), unit=self.x.unit, name=self.name)

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
    def dtype(self):
        return self.x.dtype

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_):
        self._name = name_
        for c, xyz in self._xyz.items():
            xyz.name = self._name + "_" + c

    def __add__(self, other):
        return _binary_op("__add__", self, other)

    def __iadd__(self, other):
        return _binary_op("__iadd__", self, other)

    def __sub__(self, other):
        return _binary_op("__sub__", self, other)

    def __isub__(self, other):
        return _binary_op("__isub__", self, other)

    def __mul__(self, other):
        return _binary_op("__mul__", self, other)

    def __imul__(self, other):
        return _binary_op("__imul__", self, other)

    def __truediv__(self, other):
        return _binary_op("__truediv__", self, other)

    def __itruediv__(self, other):
        return _binary_op("__itruediv__", self, other)

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        return np.reciprocal(self / other)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return -(self - other)

    def __pow__(self, number):
        return self.__class__(**{c: xyz**number for c, xyz in self._xyz.items()})

    def __neg__(self):
        return self.__class__(**{c: -xyz for c, xyz in self._xyz.items()})

    def __lt__(self, other):
        return _binary_op("__lt__", self, other)

    def __le__(self, other):
        return _binary_op("__le__", self, other)

    def __gt__(self, other):
        return _binary_op("__gt__", self, other)

    def __ge__(self, other):
        return _binary_op("__ge__", self, other)

    def __eq__(self, other):
        return _binary_op("__eq__", self, other)

    def __ne__(self, other):
        return _binary_op("__ne__", self, other)

    def __and__(self, other):
        return _binary_op("__and__", self, other)

    def __or__(self, other):
        return _binary_op("__or__", self, other)

    def __xor__(self, other):
        return _binary_op("__xor__", self, other)

    def __invert__(self):
        return np.logical_not(self)

    def to(self, unit):
        return self.__class__(**{c: xyz.to(unit) for c, xyz in self._xyz.items()})

    def _wrap_numpy(self, func, *args, **kwargs):
        if isinstance(args[0], (tuple, list)):
            # Case where we have a sequence of vectors, e.g. `concatenate`
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

    def cross(self, other):
        x = self.y * other.z
        x -= self.z * other.y
        y = self.z * other.x
        y -= self.x * other.z
        z = self.x * other.y
        z -= self.y * other.x
        return self.__class__(x, y, z)

    def _get_parent_pos(self):
        group = self.parent
        return group.get("position", group.parent["amr"]["position"])

    @property
    def r(self):
        if self.unit.is_compatible_with("meter"):
            v = self.norm
            v.name = self.name + "_r"
            return v
        pos = self._get_parent_pos()
        colatitude = _get_colatitude(pos)
        azimuth = _get_azimuth(pos)
        vec_r = np.cos(azimuth) * np.sin(colatitude) * self.x + np.sin(
            azimuth) * np.sin(colatitude) * self.y + np.cos(colatitude) * self.z
        vec_r.name = self.name + "_r"
        return vec_r

    @property
    def theta(self):
        if self.unit.is_compatible_with("meter"):
            return _get_colatitude(self)
        pos = self._get_parent_pos()
        colatitude = _get_colatitude(pos)
        azimuth = _get_azimuth(pos)
        vec_theta = np.cos(colatitude) * np.cos(azimuth) * self.x + np.cos(
            colatitude) * np.sin(azimuth) * self.y - np.sin(colatitude) * self.z
        vec_theta.name = self.name + "_theta"
        return vec_theta

    @property
    def phi(self):
        if self.unit.is_compatible_with("meter"):
            return _get_azimuth(self)
        pos = self._get_parent_pos()
        azimuth = _get_azimuth(pos)
        vec_phi = -np.sin(azimuth) * self.x + np.cos(azimuth) * self.y
        vec_phi.name = self.name + "_phi"
        return vec_phi

    @property
    def cyl_r(self):
        if self.unit.is_compatible_with("meter"):
            v = Array(values=np.linalg.norm([self.x.values, self.y.values], axis=0),
                      unit=self.unit)
            v.name = self.name + "_cyl_r"
            return v
        pos = self._get_parent_pos()
        azimuth = _get_azimuth(pos)
        vec_cyl_r = np.cos(azimuth) * self.x + np.sin(azimuth) * self.y
        vec_cyl_r.name = self.name + "_cyl_r"
        return vec_cyl_r
