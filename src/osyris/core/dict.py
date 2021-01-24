import numpy as np
from .array import Array
from ..utils import value_to_string


class Dict:
    def __init__(self):
        self._container = {}
        self.meta = {}

    def __iter__(self):
        return self._container.__iter__()

    def __getitem__(self, key):
        return self._container[key]

    def __setitem__(self, key, value):
        if isinstance(value, Array):
            value.name = key
            value.parent = self
            self._container[key] = value
        else:
            self._container[key] = Array(values=value, name=key, parent=self)

    def __delitem__(self, key):
        return self._container.__delitem__(key)

    def __repr__(self):
        return str(self)

    def __str__(self):
        # columns = ["Name", " Type", "  Group", " Unit", "  Min", "     Max"]
        columns = ["Name", " Type", "  Unit", "  Min", "     Max"]
        maxlen = {}
        for col in columns:
            maxlen[col] = 0
        print_list = {}
        for key, item in self._container.items():
            print_list[key] = {}
            print_list[key][columns[0]] = key
            print_list[key][columns[1]] = "vector" if item.ndim > 1 else "scalar"
            # print_list[key][columns[2]] = getattr(self,key).group
            print_list[key][columns[2]] = "[{}]".format(item.unit)
            print_list[key][columns[3]] = value_to_string(np.nanmin(item.values))
            print_list[key][columns[4]] = value_to_string(np.nanmax(item.values))
            for col in columns:
                maxlen[col] = max(maxlen[col], len(print_list[key][col]))

        rule = "-" * (7 + np.sum(list(maxlen.values())))
        output = ""
        # if full:
        #     for key in sorted(self.info.keys()):
        #         theShape = np.shape(self.info[key])
        #         if len(theShape) > 0:
        #             try:
        #                 output += key+": ["+str(self.info[key][0])+" ... "+str(self.info[key][-1])+"]\n"
        #             except IndexError:
        #                 output += key+": "+str(self.info[key])+"\n"
        #         else:
        #             output += key+": "+str(self.info[key])+"\n"
        #     output += "\n"
        # output += "The variables are:\n"
        for col in columns:
            output += col.ljust(maxlen[col])
        output += "\n"
        for key in sorted(print_list.keys()):
            for col in columns:
                output += print_list[key][col].ljust(maxlen[col])+" "
            output += "\n"
        return output

    def keys(self):
        return self._container.keys()

    def items(self):
        return self._container.items()

    def values(self):
        return self._container.values()