class Array:
    def __init__(self, values=None,unit=1,
        # label=None,norm=1.0,
             parent=None,name=""):

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

    @property
    def values(self):
        return self._array

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
        if self._ndim > 1:
            return self._array[:, 0]

    @property
    def y(self):
        if self._ndim > 1:
            return self._array[:, 1]

    @property
    def z(self):
        if self._ndim > 1:
            return self._array[:, 2]


    def _expect_same_unit_kind(self, other):
        if self.unit != other.unit:
            raise RuntimeError("Units don't agree in operation.")
        if self.kind != other.kind:
            raise RuntimeError("Operands are not of the same kind.")

    def _get_parent(self, other):
        if self.parent == other.parent:
            return self.parent
        else:
            return

    def _get_kind(self, other):
        return


    def __add__(self, other):
        self._expect_same_unit_kind(other)
        parent = self._get_parent(other)

        return Array(values=self.values+other.values,unit=self.unit,
             parent=parent,kind=self.kind)

    def __sub__(self, other):
        self._expect_same_unit_kind(other)
        parent = self._get_parent(other)
        return Array(values=self.values-other.values,unit=self.unit,
             parent=parent,kind=self.kind)

    def __mul__(self, other):
        if isinstance(other, Array):
            parent = self._get_parent(other)
            return Array(values=self.values*other.values,unit=self.unit * other.unit,
                 parent=parent, kind=self.kind)
        else:
            return Array(values=self.values*other, unit=self.unit, label=self.label,
            parent=self.parent, kind=self.kind, name=self.name)

    def __truediv__(self, other):
        if isinstance(other, Array):
            parent = self._get_parent(other)
            return Array(values=self.values/other.values,unit=self.unit / other.unit,
                 parent=parent, kind=self.kind)
        else:
            return Array(values=self.values/other, unit=self.unit, label=self.label,
            parent=self.parent, kind=self.kind, name=self.name)

    def __rmul__(self, number):
        return Array(values=self.values*number, unit=self.unit, label=self.label,
            parent=self.parent, kind=self.kind, name=self.name)

    def __rtruediv__(self, other):
        return Array(values=self.values/number, unit=self.unit, label=self.label,
            parent=self.parent, kind=self.kind, name=self.name)