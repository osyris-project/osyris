# SPDX-License-Identifier: BSD-3-Clause

import re


class UnitsLibrary:
    """
    A small helper class that holds the list of units, and behaves a bit like a
    dict but also performs regex matching on key with wildcard character '*'.
    """

    def __init__(self, library, default_unit):
        self._library = library
        self._default_unit = default_unit

    def __getitem__(self, key):
        if key in self._library:
            return self._library[key]
        for name, unit in self._library.items():
            if "*" in name:
                regex = re.compile(name.replace("*", ".+"))
                if re.match(regex, key):
                    return unit
        return self._default_unit

    def __setitem__(self, key, value):
        self._library[key] = value

    def __str__(self):
        return str(self._library)

    def __repr__(self):
        return str(self)

    def keys(self):
        return self._library.keys()

    def values(self):
        return self._library.values()

    def items(self):
        return self._library.items()

    def update(self, *args, **kwargs):
        self._library.update(*args, **kwargs)
