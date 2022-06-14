# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

import numpy as np


def change_origin(dataset, new_origin):
    """
    Translate all positionnal coordinates to new origin
    """
    for g in dataset.groups.keys():
        for element in dataset[g]:
            if "position" in element.lower():
                dataset[g][element] -= new_origin
    dataset.origin = new_origin


def rotation_matrix(vec, angle):
    """
    Returns 3D rotation matrix of angle 'angle' around rotation vector 'vec'.
    """
    if isinstance(vec, list):
        vec = np.array(vec)
    vec = vec / np.linalg.norm(vec)
    R = np.cos(angle) * np.identity(3) + (np.sin(angle)) * np.cross(
        vec,
        np.identity(vec.shape[0]) * -1) + (1 - np.cos(angle)) * (np.outer(vec, vec))
    return R


def get_ang_mom(subdomain, dr_L):
    """
    Compute angular momentum vector in sphere of radius dr_L
    """
    sphere = (subdomain["amr"]["position"].norm <= dr_L).values
    pos = subdomain["amr"]["position"][sphere]
    mv = subdomain["hydro"]["mass"][sphere] * subdomain["hydro"]["velocity"][sphere]
    L = np.sum(pos.cross(mv))
    return L


def change_basis(subdomain, new_basis):
    """
    Rotates all vectors in dataset to align with vector in 'new_basis'
    """
    try:
        old_basis = subdomain.basis
    except AttributeError:
        old_basis = [0, 0, 1]  # assume it's the ramses grid
    if hasattr(new_basis, 'nvec'):  # if it's a vector
        new_basis = [new_basis.x.values, new_basis.y.values, new_basis.z.values]
    rot_angle = np.arccos(
        np.dot(old_basis, new_basis) /
        (np.linalg.norm(old_basis) * np.linalg.norm(new_basis)))
    rot_vector = np.cross(new_basis, old_basis)
    R = rotation_matrix(rot_vector, rot_angle)
    for g in subdomain.groups.keys():
        for element in subdomain[g]:
            if hasattr(subdomain[g][element], 'nvec'):
                # all of this will be simplified once matmul is integraded into Vector
                vector = np.array([
                    subdomain[g][element].x.values, subdomain[g][element].y.values,
                    subdomain[g][element].z.values
                ])
                vector = R @ vector
                subdomain[g][element] = subdomain[g][element].__class__(
                    x=vector[0],
                    y=vector[1],
                    z=vector[2],
                    unit=subdomain[g][element].unit)
    subdomain.basis = new_basis
