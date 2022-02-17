# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/nvaytet/osyris)

import numpy as np
from ..core import Array
from ..config import units


def rotation_matrix(alpha, beta, gamma):
    """
    Returns the 3D rotation matrix of angles 'alpha' (yaw) around x axis
                                             'beta' (pitch) around y axis
                                             'gamma' (roll) around z axis
    """
    Rx = np.array([[1, 0, 0], [0, np.cos(alpha), -np.sin(alpha)],
                   [0, np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0],
                   [-np.sin(beta), 0, np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                   [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]])

    return Rz @ Ry @ Rx


def rotation_matrix_axis_angle(vec, angle):
    """
    Returns 3D rotation matrix of angle 'angle' around rotation vector 'vec'.
    """
    if isinstance(vec, list):
        vec = np.array(vec)
    R = np.cos(angle) * np.identity(3) + (np.sin(angle)) * np.cross(
        vec,
        np.identity(vec.shape[0]) * -1) + (1 - np.cos(angle)) * (np.outer(vec, vec))
    return R


def spherical_to_cartesian(r, theta, phi):
    """
    Converts spherical components radius, colatitude and azimuth to x,y,z cartesian coordinates
    Returns an osyris array of shape (len(r), 3) containing a stacked xyz
    """
    x, y, z = (r * np.cos(phi) * np.sin(theta), r * np.sin(phi) * np.sin(theta),
               r * np.cos(theta))
    xyz = np.vstack([x, y, z])
    return Array(values=np.transpose(xyz).values, unit=r.unit, name="position")


def cartesian_to_spherical(position, origin=[0, 0, 0]):
    """
    Converts cartesian components in 'position' to spherical components
    Returns a list of osyris Arrays containing radius, colatitude and azimuth
    """
    if isinstance(origin, Array):
        centered_pos = position - origin
    else:
        centered_pos = position - Array(values=origin, unit=centered_pos.unit)

    X = centered_pos[:, 0].values
    Y = centered_pos[:, 1].values
    Z = centered_pos[:, 2].values

    radius = np.linalg.norm([X, Y, Z], axis=0)
    colatitude = np.arctan2(np.linalg.norm([X, Y], axis=0), Z)

    azimuth = np.arctan2(Y, X) + np.pi

    radius = Array(values=radius, unit=position.unit, name='radius')
    colatitude = Array(values=colatitude, unit=units('rad'), name='colatitude')
    azimuth = Array(values=azimuth, unit=units('rad'), name='azimuth')

    return radius, colatitude, azimuth


def cartesian_to_cylindrical(position, origin=[0, 0, 0]):
    """
    Converts cartesian components in 'position' to cylindrical components
    Returns a list of osyris Arrays containing radius, azimuth and elevation
    """
    if isinstance(origin, Array):
        centered_pos = position - origin
    else:
        centered_pos = position - Array(values=origin, unit=centered_pos.unit)

    elevation = centered_pos[:, 2].values

    radius = np.linalg.norm([centered_pos[:, 0].values, centered_pos[:, 1].values],
                            axis=0)
    azimuth = np.arctan2(centered_pos[:, 1].values, centered_pos[:, 0].values)

    radius = Array(values=radius, unit=centered_pos.unit, name='radius')
    azimuth = Array(values=azimuth, unit=units('rad'), name='azimuth')
    elevation = Array(values=elevation, unit=centered_pos.unit, name='elevation')

    return radius, azimuth, elevation


def cylindrical_to_cartesian(r, azimuth, elevation):
    """
    Converts cylindrical components radius, azimuth and elevation to x,y,z cartesian coordinates
    Returns an osyris array of shape (len(r), 3) containing a stacked xyz
    """
    x, y, z = (radius * np.cos(azimuth), radius * np.sin(azimuth), elevation)
    xyz = np.vstack([x, y, z])
    return Array(values=np.transpose(xyz).values, unit=r.unit, name="position")
