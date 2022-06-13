# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)
"""
Define default values so that you don't have to specify them every time.
"""

from math import sqrt, pi


def configure_constants(units):

    units.define('bolometric_luminosity = 3.0128e+28 * W = L_bol0')
    units.define('solar_luminosity = 3.828e+26 * W = L_sun = L_sol')
    units.define('earth_mass = 5.97216787e+27 * g = M_earth')
    units.define('jupiter_mass = 1.8981246e+30 * g = M_jup')
    units.define('solar_mass = 1.9889e+33 * g = M_sun = M_sol')
    units.define('earth_radius = 6.3781e+08 * cm = R_earth')
    units.define('jupiter_radius = 7.1492e+09 * cm = R_jup')
    units.define('solar_radius = 6.957e+10 * cm = R_sun = R_sol')
    units.define(
        'radiation_constant = 7.56591469318689378e-015 * erg / cm^3 / K^4 = ar')


def configure_units(units, unit_d, unit_l, unit_t):

    density = unit_d * units("g / cm**3")
    velocity = (unit_l / unit_t) * units("cm / s")
    magnetic_field = sqrt(4.0 * pi * unit_d * (unit_l / unit_t)**2) * units("G")
    momentum = density * velocity
    acceleration = (unit_l / unit_t**2) * units("cm / s**2")
    energy = unit_d * ((unit_l / unit_t)**2) * units("erg / cm**3")
    time = unit_t * units("s")
    length = unit_l * units("cm")
    mass = density * length**3
    temperature = 1.0 * units("K")
    grav_potential = velocity**2

    library = {
        'unit_d': unit_d,
        'unit_l': unit_l,
        'unit_t': unit_t,
        'density': density,
        'velocity': velocity,
        'velocity_*': velocity,
        'momentum': momentum,
        'momentum_*': momentum,
        'magnetic_field': magnetic_field,
        'B_left': magnetic_field,
        'B_left_*': magnetic_field,
        'B_right': magnetic_field,
        'B_right_*': magnetic_field,
        'B_field': magnetic_field,
        'B_field_*': magnetic_field,
        'B_*_left': magnetic_field,
        'B_*_right': magnetic_field,
        'acceleration': acceleration,
        'grav_acceleration': acceleration,
        'grav_acceleration_*': acceleration,
        'grav_potential': grav_potential,
        'energy': energy,
        'internal_energy': energy,
        'thermal_pressure': energy,
        'pressure': energy,
        'radiative_energy': energy,
        'radiative_energy_*': energy,
        'time': time,
        'length': length,
        'x': length,
        'y': length,
        'z': length,
        'position': length,
        'position_*': length,
        'dx': length,
        'mass': mass,
        'temperature': temperature
    }
    return library


def additional_variables(data):
    """
    Here are some additional variables that are to be computed every time data
    is loaded.

    It is recommended to place your variables in a `try/except` block, which
    will prevent errors if the variables are not found, for instance when
    loading data from a different simulation.
    """

    # Magnetic field
    try:
        data['hydro']['B_field'] = 0.5 * (data['hydro']['B_left'] +
                                          data['hydro']['B_right'])
    except KeyError:
        pass

    # Mass
    try:
        data['hydro']['mass'] = (data['hydro']['density'] *
                                 data['amr']['dx']**3).to('M_sun')
    except KeyError:
        pass
