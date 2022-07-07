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


def _rotation_matrix(vec, angle):
    """
    Returns 3D rotation matrix of angle 'angle' around rotation vector 'vec'.
    """
    if isinstance(vec, list):
        vec = np.array(vec)
    vec = vec / np.linalg.norm(vec)
    r = np.cos(angle) * np.identity(3) + (np.sin(angle)) * np.cross(
        vec,
        np.identity(vec.shape[0]) * -1) + (1 - np.cos(angle)) * (np.outer(vec, vec))
    return r


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
        else:
            raise RuntimeError("Unknown basis keyword '{}'.\nAvailable keywords"
                               " are 'top' and 'side'.".format(basis))
    return basis
