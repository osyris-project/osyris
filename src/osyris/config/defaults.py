# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
"""
Define default values so that you don't have to specify them every time.
"""
import math
from pint import UnitRegistry

ureg = UnitRegistry(system="cgs")

parameters = {
    'scale': 'au',
    'path': None,
    'select': None,
    'cmap': 'viridis',
    'render_mode': 'pcolormesh'
}


def get_unit(string, ud, ul, ut):

    density = ud * (ureg.g / (ureg.cm**3))
    velocity = (ul / ut) * (ureg.cm / ureg.s)
    magnetic_field = math.sqrt(4.0 * math.pi * ud * (ul / ut)**2) * ureg.G
    momentum = density * velocity
    acceleration = (ul / ut**2) * (ureg.cm / (ureg.s**2))
    energy = ud * ((ul / ut)**2) * (ureg.erg / (ureg.cm**3))
    time = ut * ureg.s
    length = ul * ureg.cm

    ramses_units = {
        'density': density,
        'velocity': velocity,
        'velocity_x': velocity,
        'velocity_y': velocity,
        'velocity_z': velocity,
        'momentum': momentum,
        'momentum_x': momentum,
        'momentum_y': momentum,
        'momentum_z': momentum,
        'B_left': magnetic_field,
        'B_left_x': magnetic_field,
        'B_left_y': magnetic_field,
        'B_left_z': magnetic_field,
        'B_right': magnetic_field,
        'B_right_x': magnetic_field,
        'B_right_y': magnetic_field,
        'B_right_z': magnetic_field,
        'B_field': magnetic_field,
        'B_field_x': magnetic_field,
        'B_field_y': magnetic_field,
        'B_field_z': magnetic_field,
        'B_x_left': magnetic_field,
        'B_y_left': magnetic_field,
        'B_z_left': magnetic_field,
        'B_x_right': magnetic_field,
        'B_y_right': magnetic_field,
        'B_z_right': magnetic_field,
        'grav_acceleration': acceleration,
        'grav_acceleration_x': acceleration,
        'grav_acceleration_y': acceleration,
        'grav_acceleration_z': acceleration,
        'thermal_pressure': energy,
        'pressure': energy,
        'radiative_energy': energy,
        'radiative_energy_1': energy,
        'temperature': 1.0 * ureg.K,
        'time': time,
        'x': length,
        'y': length,
        'z': length,
        'dx': length
    }

    if string in ramses_units:
        return ramses_units[string]
    else:
        return 1.0 * ureg.dimensionless


def additional_units():
    """
    Define additional useful ureg and constants
    """
    ureg.define('solar_mass = 1.9889e+33 * g = msun')
    ureg.define('radiation_constant = 7.56591469318689378e-015 * erg / cm^3 / K^4 = ar')


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
        data['hydro']['mass'] = data['hydro']['density'] * data['amr']['dx']**3
        data['hydro']['mass'].to('msun')
    except KeyError:
        pass
