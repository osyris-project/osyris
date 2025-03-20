# SPDX-License-Identifier: BSD-3-Clause
import numpy as np
from pint import Quantity

from .. import units
from .array import Array
from .base import Base
from .tools import value_to_string


def perpendicular_vector(v):
    """
    Compute a vector perpendicular to the input vector
    """

    if v.z.values == 0:
        return Vector(-v.y.values, v.x.values, 0, unit=v.unit)
    else:
        return Vector(1.0, 1.0, (-1.0 * (v.x + v.y) / v.z).values, unit=v.unit)


def normalize(v):
    """
    Normalize the input vector
    """
    norm = v.norm
    nvals = norm.values
    if norm.shape:
        return v / np.where(nvals == 0, 1, nvals)
    return v / (nvals or 1)


class Vector(Base):
    def __init__(self, x, y=None, z=None, name="", unit=None):
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
        self.name = name

    def _validate_component(self, array, shape, unit):
        if array is None:
            return
        if array.shape != shape:
            raise ValueError(
                "The shape of the component does not match the "
                "shape of the x component"
            )
        if array.unit != unit:
            raise ValueError(
                "The unit of the component does not match the "
                "unit of the x component"
            )
        return array.values

    @property
    def _xyz(self):
        out = {"x": self.x}
        if self.y is not None:
            out["y"] = self.y
        if self.z is not None:
            out["z"] = self.z
        return out

    def __getitem__(self, slice_):
        return self.__class__(
            **{c: xyz[slice_] for c, xyz in self._xyz.items()}, name=self._name
        )

    def __len__(self):
        return len(self.x)

    def __str__(self):
        name_str = "'" + self._name + "' "
        xyz = self._xyz
        comps_str = ", {" + ",".join(x for x in xyz) + "}"
        if len(self) == 0:
            values_str = "Value: " + ", ".join(
                value_to_string(x.values) for x in xyz.values()
            )
            unit_str = " [{:~}] ".format(self.unit)
            shape_str = str(self.shape)
            return name_str + values_str + unit_str + shape_str + comps_str
        else:
            return str(self.norm) + comps_str

    def copy(self):
        """
        Create a (deep) copy of the vector.
        """
        return self.__class__(
            **{c: xyz.copy() for c, xyz in self._xyz.items()}, name=str(self._name)
        )

    @property
    def values(self):
        """
        A numpy array with the values of each component
        """
        return np.array([xyz.values for c, xyz in self._xyz.items()]).T

    @property
    def norm(self):
        """
        Compute the norm of the vector.
        """
        if (self.y is None) and (self.z is None):
            return self.x
        out = self.x.values * self.x.values
        out += self.y.values * self.y.values
        if self.z is not None:
            out += self.z.values * self.z.values
        return Array(values=np.sqrt(out), unit=self.x.unit, name=self.name)

    @property
    def unit(self):
        """
        The unit of the vector.
        """
        return self.x.unit

    @unit.setter
    def unit(self, unit_):
        u = units(unit_)
        self.x.unit = u
        if self.y is not None:
            self.y.unit = u
        if self.z is not None:
            self.z.unit = u

    @property
    def ndim(self):
        """
        The number of dimensions of the vector array
        (this is not the same as the number of components).
        """
        return self.x.ndim

    @property
    def nvec(self):
        """
        The number of components of the vector.
        """
        if (self.y is None) and (self.z is None):
            return 1
        if self.z is None:
            return 2
        return 3

    @property
    def shape(self):
        """
        The shape of the vector array
        (this is not the same as the shape of the components).
        """
        return self.x.shape

    @property
    def dtype(self):
        """
        The dtype of the vector array.
        """
        return self.x.dtype

    @property
    def name(self):
        """
        The name of the vector.
        """
        return self._name

    @name.setter
    def name(self, name_):
        self._name = name_
        for c, xyz in self._xyz.items():
            xyz.name = self._name + "_" + c

    def _to_vector(self, x):
        if isinstance(x, (int, float, np.ndarray, Quantity)):
            x = Array(values=x)
        if isinstance(x, Array):
            x = self.__class__(**{c: x for c in self._xyz.keys()})
        return x

    def _binary_op(self, op, other):
        other = self._to_vector(other)
        if self.nvec != other.nvec:
            raise ValueError("Operands do not have the same number of components.")
        return self.__class__(
            **{c: getattr(xyz, op)(getattr(other, c)) for c, xyz in self._xyz.items()}
        )

    def __add__(self, other):
        return self._binary_op("__add__", other)

    def __iadd__(self, other):
        return self._binary_op("__iadd__", other)

    def __sub__(self, other):
        return self._binary_op("__sub__", other)

    def __isub__(self, other):
        return self._binary_op("__isub__", other)

    def __mul__(self, other):
        return self._binary_op("__mul__", other)

    def __imul__(self, other):
        return self._binary_op("__imul__", other)

    def __truediv__(self, other):
        return self._binary_op("__truediv__", other)

    def __itruediv__(self, other):
        return self._binary_op("__itruediv__", other)

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
        return self._binary_op("__lt__", other)

    def __le__(self, other):
        return self._binary_op("__le__", other)

    def __gt__(self, other):
        return self._binary_op("__gt__", other)

    def __ge__(self, other):
        return self._binary_op("__ge__", other)

    def __eq__(self, other):
        return self._binary_op("__eq__", other)

    def __ne__(self, other):
        return self._binary_op("__ne__", other)

    def __and__(self, other):
        return self._binary_op("__and__", other)

    def __or__(self, other):
        return self._binary_op("__or__", other)

    def __xor__(self, other):
        return self._binary_op("__xor__", other)

    def __invert__(self):
        return np.logical_not(self)

    def to(self, unit):
        """
        Convert the vector to a new unit.
        """
        return self.__class__(**{c: xyz.to(unit) for c, xyz in self._xyz.items()})

    def _maybe_vector_attr(self, arg, attr):
        if isinstance(arg, (tuple, list)):
            return tuple(self._maybe_vector_attr(a, attr) for a in arg)
        if isinstance(arg, self.__class__):
            return getattr(arg, attr)
        return arg

    def _wrap_numpy(self, func, *args, **kwargs):
        if isinstance(args[0], (tuple, list)):
            arg0 = tuple(self._to_vector(arg) for arg in args[0])
        else:
            arg0 = self._to_vector(args[0])
        _args = (arg0, *[self._to_vector(arg) for arg in args[1:]])
        _kwargs = {key: self._to_vector(val) for key, val in kwargs.items()}
        return self.__class__(
            **{
                c: func(
                    *[self._maybe_vector_attr(a, c) for a in _args],
                    **{
                        k: self._maybe_vector_attr(kwa, c) for k, kwa in _kwargs.items()
                    },
                )
                for c in self._xyz.keys()
            }
        )

    def reshape(self, *shape):
        """
        Reshape the vector arrays.

        Parameters
        ----------
        shape : tuple
            The new shape of the vector arrays.
        """
        return self.__class__(
            **{c: xyz.reshape(*shape) for c, xyz in self._xyz.items()}
        )

    @property
    def nbytes(self):
        """
        The number of bytes used by the vector.
        """
        return np.sum([xyz.nbytes for xyz in self._xyz.values()])

    def dot(self, other):
        """
        Compute the dot product of two vectors.
        """
        out = np.zeros(self.shape)
        for c1, c2 in zip(self._xyz.values(), other._xyz.values()):
            out += (c1 * c2).values
        return Array(values=out, unit=self.unit * other.unit)

    def cross(self, other):
        """
        Compute the cross product of two vectors.
        """
        x = self.y * other.z
        x -= self.z * other.y
        y = self.z * other.x
        y -= self.x * other.z
        z = self.x * other.y
        z -= self.y * other.x
        return self.__class__(x, y, z)


class VectorBasis:
    def __init__(self, n, u=None, v=None):
        self.n = n
        self.u = perpendicular_vector(self.n) if u is None else u
        self.v = self.n.cross(self.u) if v is None else v
        self.n = normalize(self.n)
        self.u = normalize(self.u)
        self.v = normalize(self.v)
        self.n.name = n.name
        if u is not None:
            self.u.name = u.name
        if v is not None:
            self.v.name = v.name

    def roll(self):
        return VectorBasis(n=self.u, u=self.v, v=self.n)

    def __str__(self) -> str:
        header = "VectorBasis:\n    "
        body = "\n    ".join(f"{x}: {getattr(self, x)}" for x in ("n", "u", "v"))
        return header + body

    def __repr__(self) -> str:
        return str(self)
