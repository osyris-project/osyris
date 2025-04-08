# SPDX-License-Identifier: BSD-3-Clause

import re
import numpy as np

from .units import units


class UnitsLibrary:
    """
    A small helper class that holds the list of units, and behaves a bit like a
    dict but also performs regex matching on key with wildcard character '*'.
    """

    def __init__(self, *, unit_d, unit_l, unit_t, default_unit):
        density = unit_d * units("g / cm**3")
        velocity = (unit_l / unit_t) * units("cm / s")
        magnetic_field = np.sqrt(4.0 * np.pi * unit_d * (unit_l / unit_t) ** 2) * units(
            "G"
        )
        momentum = density * velocity
        acceleration = (unit_l / unit_t**2) * units("cm / s**2")
        energy = unit_d * ((unit_l / unit_t) ** 2) * units("erg / cm**3")
        time = unit_t * units("s")
        length = unit_l * units("cm")
        mass = density * length**3
        temperature = 1.0 * units("K")
        grav_potential = velocity**2

        self._library = {
            "unit_d": unit_d,
            "unit_l": unit_l,
            "unit_t": unit_t,
            "density": density,
            "velocity": velocity,
            "velocity_*": velocity,
            "momentum": momentum,
            "momentum_*": momentum,
            "magnetic_field": magnetic_field,
            "B_left": magnetic_field,
            "B_left_*": magnetic_field,
            "B_right": magnetic_field,
            "B_right_*": magnetic_field,
            "B_field": magnetic_field,
            "B_field_*": magnetic_field,
            "B_*_left": magnetic_field,
            "B_*_right": magnetic_field,
            "acceleration": acceleration,
            "grav_acceleration": acceleration,
            "grav_acceleration_*": acceleration,
            "grav_potential": grav_potential,
            "energy": energy,
            "internal_energy": energy,
            "thermal_pressure": energy,
            "pressure": energy,
            "radiative_energy": energy,
            "radiative_energy_*": energy,
            "time": time,
            "length": length,
            "x": length,
            "y": length,
            "z": length,
            "position": length,
            "position_*": length,
            "dx": length,
            "mass": mass,
            "temperature": temperature,
        }
        # self._library = library
        self._default_unit = default_unit

    def __getitem__(self, key):
        if key in self._library:
            return self._library[key]
        for name, unit in self._library.items():
            if "*" in name:
                regex = re.compile(name.replace("*", ".+"))
                if re.match(regex, key):
                    return unit
        return self._default_unit

    def __setitem__(self, key, value):
        self._library[key] = value

    def __str__(self):
        return str(self._library)

    def __repr__(self):
        return str(self)

    def keys(self):
        return self._library.keys()

    def values(self):
        return self._library.values()

    def items(self):
        return self._library.items()

    def update(self, *args, **kwargs):
        self._library.update(*args, **kwargs)
