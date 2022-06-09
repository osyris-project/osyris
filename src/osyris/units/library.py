# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import re
from pint import Unit


class UnitsLibrary:
    """
    A small helper class that holds the list of units, and behaves a bit like a
    dict but also performs regex matching on key with wildcard character '*'.
    """
    def __init__(self, *args, **kwargs):
        self._library = dict(*args, **kwargs)

    def __getitem__(self, key):
        if key in self._library:
            return self._library[key]
        for name, unit in self._library.items():
            if '*' in name:
                regex = re.compile(name.replace('*', '.+'))
                if re.match(regex, key):
                    return unit
        return 1.0 * Unit('')

    def __setitem__(self, key, value):
        self._library[key] = value

    def update(self, *args, **kwargs):
        self._library.update(*args, **kwargs)
