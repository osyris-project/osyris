# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

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


def cartesian_to_spherical_rotation_matrix(colatitude, azimuth):
    """
    Used to convert arbitrary cartesian vectors into spherical vectors
    """

    R = np.array([[
        np.sin(colatitude) * np.cos(azimuth),
        np.sin(colatitude) * np.sin(azimuth),
        np.cos(colatitude)
    ],
                  [
                      np.cos(colatitude) * np.cos(azimuth),
                      np.cos(colatitude) * np.sin(azimuth), -np.sin(colatitude)
                  ], [-np.sin(azimuth), np.cos(azimuth), 0]],
                 dtype=object)
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


def get_spherical_colatitude(position_vec):
    """
    Compute colatitude from cartesian position vector
    """
    X = position_vec.x.values
    Y = position_vec.y.values
    Z = position_vec.z.values

    theta = np.arctan2(np.linalg.norm([X, Y], axis=0), Z)

    return position_vec.x.__class__(values=theta, unit='rad', name="position_phi")


def get_spherical_azimuth(position_vec):
    """
    Compute colatitude from cartesian position vector
    """
    X = position_vec.x.values
    Y = position_vec.y.values
    phi = np.arctan2(Y, X)

    return position_vec.x.__class__(values=phi, unit='rad', name="position_phi")


def get_spherical_components(vec, comp):
    """
    Get the spherical components of an arbitrary cartesian vector 'vec'
    """
    if not hasattr(vec, 'nvec'):
        raise ValueError("This is not a vector")
    if vec.parent.name == "hydro":
        position = vec.parent.parent["amr"]["position"]
    else:
        position = vec.parent["position"]

    theta = position.theta
    phi = position.phi

    if comp.lower() == "radius":
        vec_r = np.cos(phi) * np.sin(theta) * vec.x + np.sin(phi) * np.sin(
            theta) * vec.y + np.cos(theta) * vec.z
        vec_r.name = vec.name + "_r"
        return vec_r
    elif comp.lower() == "colatitude":
        vec_theta = np.cos(theta) * np.cos(phi) * vec.x + np.cos(theta) * np.sin(
            phi) * vec.y - np.sin(theta) * vec.z
        vec_theta.name = vec.name + "_theta"
        return vec_theta
    elif comp.lower() == "azimuth":
        vec_phi = -np.sin(phi) * vec.x + np.cos(phi) * vec.y
        vec_phi.name = vec.name + "_phi"
        return vec_phi
    else:
        raise ValueError(
            "Urecognized component keyword '{}'. Valid strings are 'radius', 'colatitude' and 'azimuth'"
        )
