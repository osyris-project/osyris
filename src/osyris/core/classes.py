
class Plot:

    def __init__(self, x, y, layers, fig, ax):
        self.x = x
        self.y = y
        self.layers = layers
        self.fig = fig
        self.ax = ax


class Array:
    def __init__(self, values=None,unit=1,label=None,norm=1.0,
             parent=None,kind="scalar",name=""):

        self._values = values
        self._unit = unit
        # self._label = label
        # self.operation = operation
        # self.depth = depth
        # self.norm = norm
        self.kind = "vector" if values.ndim > 1 else "scalar"
        self.parent = parent
        self.name = name
        # self.group = group
        # self.vector = vector

    @property
    def foo(self):
        return self._foo

    @property
    def foo(self):
        return self._foo


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
        parent = self._get_parent(other)
        return Array(values=self.values/other.values,unit=self.unit / other.unit,
             parent=parent, kind=self.kind)

    def __rmul__(self, number):
        return Array(values=self.values*number, unit=self.unit, label=self.label,
            parent=self.parent, kind=self.kind, name=self.name)


class Dict:
    def __init__(self):
        self.container = {}
        self.meta = {}

    def __iter__(self):
        return self.container.__iter__()

    def __getitem__(self, key):
        return self.container[key]

    def __setitem__(self, key, value):
        if isinstance(value, Array):
            value.name = key
            value.parent = self
            self.container[key] = value
        else:
            self.container[key] = Array(values=value, name=key, parent=self)

    def keys(self):
        return self.container.keys()

    def items(self):
        return self.container.items()

    def values(self):
        return self.container.values()



def array(**kwargs):
    return Array(**kwargs)
