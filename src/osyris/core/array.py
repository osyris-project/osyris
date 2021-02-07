import numpy as np
from pint.quantity import Quantity
from .tools import value_to_string
from .. import units

class Array:
    def __init__(self, values=None, unit=None,
        # label=None,norm=1.0,
             parent=None, name=""):

        self._array = values
        self._unit = unit if unit is not None else 1.0 * units.dimensionless
        # self._label = label
        # self.operation = operation
        # self.depth = depth
        # self.norm = norm
        # self._kind = "vector" if values.ndim > 1 else "scalar"
        # self._vector = values.ndim > 1
        self._parent = parent
        self._name = name
        # self.group = group
        # self.vector = vector

    def __array__(self):
        return self._array

    def __getitem__(self, slice_):
        return self.__class__(values=self._array[slice_], unit=self._unit,
                parent=self._parent, name=self._name)

    def __str__(self):
        name_str = "'"+self._name + "' "
        values_str = "Min: " + value_to_string(
            self.values.min()) + " Max: " + value_to_string(
            self.values.max())
        unit_str = " [{:~}] ".format(self._unit.units)
        shape_str = str(self._array.shape)
        return name_str + values_str + unit_str + shape_str

    def __repr__(self):
        return str(self)

    @property
    def values(self):
        if self._array.ndim < 2:
            return self._array
        else:
            return np.linalg.norm(self._array, axis=1)

    @values.setter
    def values(self, values_):
        self._array = values_

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, unit_):
        self._unit = unit_

    @property
    def ndim(self):
        return self._array.ndim

    @property
    def shape(self):
        return self._array.shape

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
    def x(self):
        if self.ndim > 1:
            return self.__class__(values=self._array[:, 0], unit=self._unit,
                parent=self._parent, name=self._name+"_x")

    @property
    def y(self):
        if self.ndim > 1:
            return self.__class__(values=self._array[:, 1], unit=self._unit,
                parent=self._parent, name=self._name+"_y")

    @property
    def z(self):
        if self.ndim > 1:
            return self.__class__(values=self._array[:, 2], unit=self._unit,
                parent=self._parent, name=self._name+"_z")


    # def _expect_same_unit(self, other):
    #     if self._unit != other.unit:
    #         raise RuntimeError("Units don't agree in operation.")
    #     # if self.kind != other.kind:
    #     #     raise RuntimeError("Operands are not of the same kind.")

    # def _get_parent(self, other):
    #     if self._parent == other.parent:
    #         return self._parent
    #     else:
    #         return

    @property
    def label(self, name=True, unit=True):
        lab = ""
        if name:
            lab += self._name + " "
        if unit:
            lab += "[{:~}]".format(self._unit.units)
        return lab.strip()

    # def _get_kind(self, other):
    #     return

    def _broadcast(self, lhs, rhs):
        if lhs.ndim == rhs.ndim:
            return lhs, rhs
        if lhs.ndim > rhs.ndim:
            ind = np.argmax(np.array(lhs.shape) == rhs.shape[0])
            if ind == 0:
                return lhs, rhs.reshape(rhs.shape + tuple([1]))
            else:
                return lhs, rhs.reshape(tuple([1]) + rhs.shape)
        else:
            ind = np.argmax(np.array(rhs.shape) == lhs.shape[0])
            if ind == 0:
                return lhs.reshape(lhs.shape + tuple([1])), rhs
            else:
                return lhs.reshape(tuple([1]) + lhs.shape), rhs


    def __add__(self, other):
        if isinstance(other, self.__class__):
            scale_r = other.unit.to(self._unit.units)
            lhs = self._array
            rhs = other._array*scale_r.magnitude
            lhs, rhs = self._broadcast(lhs, rhs)
            return self.__class__(values=lhs+rhs, unit=self._unit)
        if isinstance(other, Quantity):
            return self.__class__(values=self._array+other.to(self._unit.units).magnitude,
                unit=self._unit)
        raise TypeError("Could not add types {} and {}.".format(type(self), type(other)))


    def __sub__(self, other):
        if isinstance(other, self.__class__):
            scale_r = other.unit.to(self._unit.units)
            lhs = self._array
            rhs = other._array*scale_r.magnitude
            lhs, rhs = self._broadcast(lhs, rhs)
            return self.__class__(values=lhs - rhs, unit=self._unit)
        if isinstance(other, Quantity):
            return self.__class__(values=self._array-other.to(self._unit.units),
                unit=self._unit)
        raise TypeError("Could not subtract types {} and {}.".format(type(self), type(other)))


    def __mul__(self, other):
        if isinstance(other, self.__class__):
            # parent = self._get_parent(other)
            scale_l = self._unit.to_base_units()
            scale_r = other._unit.to_base_units()
            result = scale_l * scale_r
            # return self.__class__(values=self._array*scale_l.magnitude *
            #     other._array * scale_r.magnitude ,unit=1.0 * result.units,
            #      parent=parent)

            lhs = self._array
            rhs = other._array*result.magnitude
            lhs, rhs = self._broadcast(lhs, rhs)
            return self.__class__(values=lhs * rhs, unit=1.0 * result.units)
        if isinstance(other, Quantity):
            scale_l = self._unit.to_base_units()
            scale_r = other.to_base_units()
            result = scale_l * scale_r
            return self.__class__(values=self._array * result.magnitude,
                unit=1.0 * result.units)

        return self.__class__(values=self._array*other, unit=self._unit)

    def __truediv__(self, other):
        if isinstance(other, self.__class__):
            # parent = self._get_parent(other)
            scale_l = self._unit.to_base_units()
            scale_r = other._unit.to_base_units()
            result = scale_l / scale_r
            # return self.__class__(values=self._array*scale_l.magnitude *
            #     other._array * scale_r.magnitude ,unit=1.0 * result.units,
            #      parent=parent)

            lhs = self._array
            rhs = other._array / result.magnitude
            lhs, rhs = self._broadcast(lhs, rhs)

            return self.__class__(values=lhs / rhs, unit=1.0 * result.units)
        if isinstance(other, Quantity):
            scale_l = self._unit.to_base_units()
            scale_r = other.to_base_units()
            result = scale_l / scale_r
            return self.__class__(values=self._array * result.magnitude,
                unit=1.0 * result.units)

        return self.__class__(values=self._array/other, unit=self._unit)




    # def __truediv__(self, other):
    #     if isinstance(other, Array):
    #         parent = self._get_parent(other)
    #         return self.__class__(values=self._array/other._array,unit=self._unit / other.unit,
    #              parent=parent)
    #     else:
    #         return self.__class__(values=self._array/other, unit=self._unit,
    #         parent=self.parent, name=self._name)

    def __rmul__(self, other):
        return self * other
        #     # scale_l = self._unit.to_base_units()
        #     # scale_r = other.to_base_units()
        #     # result = scale_l * scale_r
        #     # return self.__class__(values=self._array * result.magnitude,
        #     #     unit=1.0 * result.units)

        # return self.__class__(values=self._array*other, unit=self._unit,
        #     name=self._name)



        # return self.__class__(values=self._array*number, unit=self._unit,
        #     parent=self.parent, name=self._name)

    # def __rtruediv__(self, other):
    #     return self.__class__(values=self._array/number, unit=self._unit,
    #         parent=self.parent, name=self._name)

    def __rtruediv__(self, other):
        if isinstance(other, self.__class__):
            # parent = self._get_parent(other)
            scale_r = self._unit.to_base_units()
            scale_l = other._unit.to_base_units()
            result = scale_l / scale_r
            # return self.__class__(values=self._array*scale_l.magnitude *
            #     other._array * scale_r.magnitude ,unit=1.0 * result.units,
            #      parent=parent)

            lhs = self._array
            rhs = other._array / result.magnitude
            lhs, rhs = self._broadcast(lhs, rhs)

            return self.__class__(values=lhs / rhs, unit=1.0 * result.units)
        if isinstance(other, Quantity):
            scale_r = self._unit.to_base_units()
            scale_l = other.to_base_units()
            result = scale_l / scale_r
            return self.__class__(values=self._array * result.magnitude,
                unit=1.0 * result.units)

        return self.__class__(values=other / self._array, unit=1.0 / self._unit)



    def to(self, unit):
        new_unit = units(unit)
        ratio = self._unit.to(new_unit) / new_unit
        self._unit = 1.0 * new_unit
        self._array *= ratio.magnitude
