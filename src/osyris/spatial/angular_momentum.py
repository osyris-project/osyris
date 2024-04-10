# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2024 Osyris contributors (https://github.com/osyris-project/osyris)

from typing import Optional
from pint import Quantity
import numpy as np

from ..core import Array, Dataset, Vector, VectorBasis
from ..core.vector import perpendicular_vector


def angular_momentum_vector(
    position: Array,
    mass: Array,
    velocity: Array,
    radius: Quantity,
    origin: Optional[Vector] = None,
) -> Vector:
    """
    Calculate the angular momentum vector of the gas contained in a sphere of radius
    ``radius`` centered around the ``origin``.

    :param position: The position of the gas particles.

    :param mass: The mass of the gas particles.

    :param velocity: The velocity of the gas particles.

    :param radius: The radius of the sphere.

    :param origin: The origin of the sphere.
    """
    if origin is not None:
        position = position - origin
    sphere = (position.norm <= radius).values
    weighted_pos = position[sphere] * mass[sphere]
    return np.sum(weighted_pos.cross(velocity[sphere]))


def top_view(
    dataset: Dataset, radius: Quantity, origin: Optional[Vector] = None
) -> VectorBasis:
    """
    Find the direction vectors for a top view of the dataset according to the local
    angular momentum.

    :param dataset: The dataset.

    :param radius: The radius of the sphere.

    :param origin: The origin of the sphere.
    """
    am = angular_momentum_vector(
        position=dataset["amr"]["position"],
        mass=dataset["hydro"]["mass"],
        velocity=dataset["hydro"]["velocity"],
        radius=radius,
        origin=origin,
    )
    am = am / am.norm
    # VectorBasis will compute other components if not supplied
    return VectorBasis(n=am)


def side_view(
    dataset: Dataset, radius: Quantity, origin: Optional[Vector] = None
) -> VectorBasis:
    """
    Find the direction vectors for a side view of the dataset according to the local
    angular momentum.

    :param dataset: The dataset.

    :param radius: The radius of the sphere.

    :param origin: The origin of the sphere.
    """
    am = angular_momentum_vector(
        position=dataset["amr"]["position"],
        mass=dataset["hydro"]["mass"],
        velocity=dataset["hydro"]["velocity"],
        radius=radius,
        origin=origin,
    )
    am = am / am.norm
    return VectorBasis(n=perpendicular_vector(am), u=am)
