# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np


def get_ang_mom(subdomain, dr_L):
    """
    Compute angular momentum vector in sphere of radius dr_L
    """
    sphere = (subdomain["amr"]["position"].norm <= dr_L).values
    pos = subdomain["amr"]["position"][sphere]
    mv = subdomain["hydro"]["mass"][sphere] * subdomain["hydro"]["velocity"][sphere]
    L = np.sum(pos.cross(mv))
    return L


def _parse_basis(subdomain, basis, dr_L):
    if isinstance(basis, str):
        if dr_L is None:
            raise ValueError("Please provide the radius size with " +
                             "which to compute angular momentum (dr_L)")
        ang_mom = get_ang_mom(subdomain, dr_L)
        if basis.lower() == "top":
            basis = ang_mom / ang_mom.norm
        elif basis.lower() == "side":
            perp_v = ang_mom.__class__(1.0,
                                       1.0, (-1.0 * (ang_mom.x + ang_mom.y) /
                                             ang_mom.z).values,
                                       unit=ang_mom.unit)
            basis = perp_v / perp_v.norm
    return basis
