# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np
from .. import Dataset
import warnings


def extract_sphere(dataset, radius, origin):
    """
    Extract a spherical subdomain around an origin point.
    """
    subdomain = Dataset()
    subdomain.meta = dataset.meta.copy()

    for name, group in dataset.items():
        pos = group.get("position", group.parent["amr"]["position"])
        if pos.shape != group.shape:
            warnings.warn(
                "Ignoring datagroup '{}', which has no position ".format(group)
                + "vector and has different shape than 'amr' group."
            )
            continue
        r = (pos - origin).norm
        c = (r < radius).values
        if np.any(c):
            subdomain[name] = dataset[name][c]

    return subdomain


def extract_box(dataset, dx, dy, dz, origin):
    """
    Extract a cubic domain of size dx, dy & dz around an origin point
    """
    subdomain = Dataset()
    subdomain.meta = dataset.meta.copy()

    for name, group in dataset.items():
        pos = group.get("position", group.parent["amr"]["position"])
        if pos.shape != group.shape:
            warnings.warn(
                "Ignoring datagroup '{}', which has no position ".format(group)
                + "vector and has different shape than 'amr' group."
            )
            continue
        centered_pos = pos - origin
        cx = (centered_pos.x <= dx * 0.5) & (centered_pos.x >= -dx * 0.5)
        cy = (centered_pos.y <= dy * 0.5) & (centered_pos.y >= -dy * 0.5)
        cz = (centered_pos.z <= dz * 0.5) & (centered_pos.z >= -dz * 0.5)
        c = (cx & cy & cz).values
        if np.any(c):
            subdomain[name] = dataset[name][c]

    return subdomain
