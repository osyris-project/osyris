# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np


def extract_sphere(dataset, radius, origin):
    """
    Extract a spherical subdomain around an origin point.
    """
    subdomain = dataset.__class__(nout=dataset.meta["nout"], path=dataset.meta["path"])
    subdomain.meta = dataset.meta
    subdomain._parent = dataset

    for name, group in dataset.items():
        pos = group.get("position", group.parent["amr"]["position"])
        r = (pos - origin).norm
        c = (r < radius).values
        if np.any(c):
            subdomain[name] = dataset[name][c]

    return subdomain


def extract_cube(dataset, dx, dy, dz, origin):
    """
    Extract a cubic domain of size dx, dy & dz around an origin point
    """
    subdomain = dataset.__class__(nout=dataset.meta["nout"], path=dataset.meta["path"])
    subdomain.meta = dataset.meta
    subdomain._parent = dataset

    for name, group in dataset.items():
        pos = group.get("position", group.parent["amr"]["position"])
        centered_pos = pos - origin
        cx = (centered_pos.x <= dx * .5) & (centered_pos.x >= -dx * .5)
        cy = (centered_pos.y <= dy * .5) & (centered_pos.y >= -dy * .5)
        cz = (centered_pos.z <= dz * .5) & (centered_pos.z >= -dz * .5)
        c = (cx & cy & cz).values
        if np.any(c):
            subdomain[name] = dataset[name][c]

    return subdomain
