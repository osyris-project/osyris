# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np


def arrayclose(a, b):
    if a.unit != b.unit:
        return False
    return np.allclose(a.values, b.values)


def arraytrue(a):
    return all(a.values)


def arrayequal(a, b):
    if a.unit != b.unit:
        return False
    return np.array_equal(a.values, b.values)


def vectorclose(a, b):
    ok = arrayclose(a.x, b.x)
    if a.y is not None:
        ok *= arrayclose(a.y, b.y)
    if a.z is not None:
        ok *= arrayclose(a.z, b.z)
    return bool(ok)


def vectortrue(a):
    ok = arraytrue(a.x)
    if a.y is not None:
        ok *= arraytrue(a.y)
    if a.z is not None:
        ok *= arraytrue(a.z)
    return bool(ok)


def vectorequal(a, b):
    ok = arrayequal(a.x, b.x)
    if a.y is not None:
        ok *= arrayequal(a.y, b.y)
    if a.z is not None:
        ok *= arrayequal(a.z, b.z)
    return bool(ok)
