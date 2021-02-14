import numpy as np
from pint.quantity import Quantity
from pint.unit import Unit
from .tools import value_to_string, make_label
from .. import units



def _comparison_operator(lhs, rhs, op):
    func = getattr(np, op)
    if isinstance(rhs, Array):
        scale_r = rhs.unit.to(lhs._unit.units)
        return func(lhs._array, rhs._array*scale_r.magnitude)
    if isinstance(rhs, Quantity):
        return func(lhs._array < rhs.to(lhs._unit.units).magnitude)
    return func(lhs._array, rhs)



class Array:
    def __init__(self, values=None, unit=None,
        # label=None,norm=1.0,
             parent=None, name=""):

        self._array = values

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
        # self._label = label
        # self.operation = operation
        # self.depth = depth
        # self.norm = norm
        # self._kind = "vector" if values.ndim > 1 else "scalar"
        # self._vector = values.ndim > 1
        self._parent = parent
        self._name = name
        if self._array.shape:
            self._ndim = self._array.shape[-1]
        else:
        	self._ndim = 0
        # self.group = group
        # self.vector = vector
        self.special_functions = {"sqrt": self._sqrt}

    # def __array__(self):
    #     return self._array

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
        return self._array

    @values.setter
    def values(self, values_):
        self._array = values_

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, array_):
        self._array = array_

    @property
    def norm(self):
        if self._array.ndim < 2:
            return self._array
        else:
            return np.linalg.norm(self._array, axis=1)

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, unit_):
        self._unit = unit_

    @property
    def ndim(self):
        return self._ndim

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
    def label(self):
        return make_label(name=self._name, unit=self._unit.units)

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

    def _raise_incompatible_units_error(other, op):
        raise TypeError("Could not {} types {} and {}.".format(
            op, type(self), type(other)))


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
        self._raise_incompatible_units_error(other, "add")

    def __iadd__(self, other):
        if isinstance(other, self.__class__):
            scale_r = other.unit.to(self._unit.units)
            # lhs = self._array
            rhs = other._array*scale_r.magnitude
            # lhs, rhs = self._broadcast(lhs, rhs)
            # return self.__class__(values=lhs+rhs, unit=self._unit)
            self._array += rhs
        elif isinstance(other, Quantity):
            self._array += other.to(self._unit.units).magnitude
            # return self.__class__(values=self._array+other.to(self._unit.units).magnitude,
            #     unit=self._unit)
        else:
            self._raise_incompatible_units_error(other, "add")
        return self


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
        self._raise_incompatible_units_error(other, "subtract")


    def __isub__(self, other):
        if isinstance(other, self.__class__):
            scale_r = other.unit.to(self._unit.units)
            rhs = other._array*scale_r.magnitude
            self._array -= rhs
        elif isinstance(other, Quantity):
            self._array -= other.to(self._unit.units).magnitude
        else:
            self._raise_incompatible_units_error(other, "subtract")
        return self


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

    def __imul__(self, other):
        if isinstance(other, self.__class__):
            scale_l = self._unit.to_base_units()
            scale_r = other._unit.to_base_units()
            result = scale_l * scale_r
            # lhs = self._array
            rhs = other._array*result.magnitude
            # lhs, rhs = self._broadcast(lhs, rhs)
            # return self.__class__(values=lhs * rhs, unit=1.0 * result.units)
            self._array *= rhs
            self._unit = 1.0 * result.units
        elif isinstance(other, Quantity):
            scale_l = self._unit.to_base_units()
            scale_r = other.to_base_units()
            result = scale_l * scale_r
            self._array *= result.magnitude
            self._unit = 1.0 * result.units
            # return self.__class__(values=self._array * result.magnitude,
            #     unit=1.0 * result.units)
        else:
            self._array *= other
        return self


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


    def __itruediv__(self, other):
        if isinstance(other, self.__class__):
            # parent = self._get_parent(other)
            scale_l = self._unit.to_base_units()
            scale_r = other._unit.to_base_units()
            result = scale_l / scale_r
            # return self.__class__(values=self._array*scale_l.magnitude *
            #     other._array * scale_r.magnitude ,unit=1.0 * result.units,
            #      parent=parent)
            # lhs = self._array
            rhs = other._array / result.magnitude
            # lhs, rhs = self._broadcast(lhs, rhs)
            self._array /= rhs
            self._unit = 1.0 * result.units
            # return self.__class__(values=lhs / rhs, unit=1.0 * result.units)
        elif isinstance(other, Quantity):
            scale_l = self._unit.to_base_units()
            scale_r = other.to_base_units()
            result = scale_l / scale_r
            self._array *= result.magnitude
            self._unit = 1.0 * result.units
            # return self.__class__(values=self._array * result.magnitude,
            #     unit=1.0 * result.units)
        else:
            self._array /= other
        return self


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

    def __lt__(self, other):
        return _comparison_operator(self, other, "less")

    def __le__(self, other):
        return _comparison_operator(self, other, "less_equal")

    def __gt__(self, other):
        return _comparison_operator(self, other, "greater")

    def __ge__(self, other):
        return _comparison_operator(self, other, "greater_equal")

    def __eq__(self, other):
        return _comparison_operator(self, other, "equal")

    def __ne__(self, other):
        return _comparison_operator(self, other, "not_equal")


    def to(self, unit):
        new_unit = units(unit)
        ratio = self._unit.to(new_unit) / new_unit
        self._unit = 1.0 * new_unit
        self._array *= ratio.magnitude


    def _wrap_numpy(self, func, *args, **kwargs):
        if func.__name__ in self.special_functions:
            unit = self.special_functions[func.__name__]()
        else:
            unit = self._unit
        args = (args[0]._array,) + args[1:]
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
        Numpy array_function protocol to allow Array to work with numpy functions.
        """
        return self._wrap_numpy(func, *args, **kwargs)


    def _sqrt(self):
        """
        Handle unit in sqrt operation.
        """
        return np.sqrt(self._unit)
