import numpy as np
from ..utils import value_to_string

class Array:
    def __init__(self, values=None, unit=None,
        # label=None,norm=1.0,
             parent=None, name=""):

        self._array = values
        self._unit = unit
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


    def __getitem__(self, slice_):
        return Array(values=self._array[slice_], unit=self._unit,
                parent=self._parent, name="")

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
            return Array(values=self._array[:, 0], unit=self._unit,
                parent=self._parent, name=self.name+"_x")

    @property
    def y(self):
        if self.ndim > 1:
            return Array(values=self._array[:, 1], unit=self._unit,
                parent=self._parent, name=self.name+"_y")

    @property
    def z(self):
        if self.ndim > 2:
            return Array(values=self._array[:, 2], unit=self._unit,
                parent=self._parent, name=self.name+"_z")


    def _expect_same_unit(self, other):
        if self.unit != other.unit:
            raise RuntimeError("Units don't agree in operation.")
        # if self.kind != other.kind:
        #     raise RuntimeError("Operands are not of the same kind.")

    def _get_parent(self, other):
        if self.parent == other.parent:
            return self.parent
        else:
            return

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


    def __add__(self, other):
        self._expect_same_unit(other)
        parent = self._get_parent(other)

        return Array(values=self._array+other._array,unit=self.unit,
             parent=parent)

    def __sub__(self, other):
        self._expect_same_unit(other)
        parent = self._get_parent(other)
        return Array(values=self._array-other._array,unit=self.unit,
             parent=parent)

    def __mul__(self, other):
        if isinstance(other, Array):
            parent = self._get_parent(other)
            return Array(values=self._array*other._array,unit=self.unit * other.unit,
                 parent=parent)
        else:
            return Array(values=self._array*other, unit=self.unit,
            parent=self.parent, name=self.name)

    def __truediv__(self, other):
        if isinstance(other, Array):
            parent = self._get_parent(other)
            return Array(values=self._array/other._array,unit=self.unit / other.unit,
                 parent=parent)
        else:
            return Array(values=self._array/other, unit=self.unit,
            parent=self.parent, name=self.name)

    def __rmul__(self, number):
        return Array(values=self._array*number, unit=self.unit,
            parent=self.parent, name=self.name)

    def __rtruediv__(self, other):
        return Array(values=self._array/number, unit=self.unit,
            parent=self.parent, name=self.name)
