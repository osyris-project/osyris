# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

from .array import Array

import numpy as np


def dot(v1, v2):
    out = np.zeros(v1.x.shape)
    for (c1, c2) in zip(v1._xyz.values(), v2._xyz.values()):
        out += (c1 * c2).values
    return Array(values=out, unit=v1.x.unit * v2.x.unit)
