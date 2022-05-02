# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np


def allclose(a, b):
    if a.unit != b.unit:
        return False
    return np.allclose(a.values, b.values)


def alltrue(a):
    return all(a.values)
