# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2022 Osyris contributors (https://github.com/osyris-project/osyris)

from pint import UnitRegistry
from math import sqrt, pi
# from . import config


class Units:
    def __init__(self):
        self._ureg = UnitRegistry(system="cgs")
        self._base_units = {}

        self._ureg.define('bolometric_luminosity = 3.0128e+28 * W = L_bol0')
        self._ureg.define('solar_luminosity = 3.828e+26 * W = L_sun = L_sol')
        self._ureg.define('earth_mass = 5.97216787e+27 * g = M_earth')
        self._ureg.define('jupiter_mass = 1.8981246e+30 * g = M_jup')
        self._ureg.define('solar_mass = 1.9889e+33 * g = M_sun = M_sol')
        self._ureg.define('earth_radius = 6.3781e+08 * cm = R_earth')
        self._ureg.define('jupiter_radius = 7.1492e+09 * cm = R_jup')
        self._ureg.define('solar_radius = 6.957e+10 * cm = R_sun = R_sol')
        self._ureg.define(
            'radiation_constant = 7.56591469318689378e-015 * erg / cm^3 / K^4 = ar')

        self._ramses_units = {
            'density': 'density',
            'velocity': 'velocity',
            'velocity_x': 'velocity',
            'velocity_y': 'velocity',
            'velocity_z': 'velocity',
            'momentum': 'momentum',
            'momentum_x': 'momentum',
            'momentum_y': 'momentum',
            'momentum_z': 'momentum',
            'B_left': 'magnetic_field',
            'B_left_x': 'magnetic_field',
            'B_left_y': 'magnetic_field',
            'B_left_z': 'magnetic_field',
            'B_right': 'magnetic_field',
            'B_right_x': 'magnetic_field',
            'B_right_y': 'magnetic_field',
            'B_right_z': 'magnetic_field',
            'B_field': 'magnetic_field',
            'B_field_x': 'magnetic_field',
            'B_field_y': 'magnetic_field',
            'B_field_z': 'magnetic_field',
            'B_x_left': 'magnetic_field',
            'B_y_left': 'magnetic_field',
            'B_z_left': 'magnetic_field',
            'B_x_right': 'magnetic_field',
            'B_y_right': 'magnetic_field',
            'B_z_right': 'magnetic_field',
            'grav_acceleration': 'acceleration',
            'grav_acceleration_x': 'acceleration',
            'grav_acceleration_y': 'acceleration',
            'grav_acceleration_z': 'acceleration',
            'thermal_pressure': 'energy',
            'pressure': 'energy',
            'radiative_energy': 'energy',
            'radiative_energy_1': 'energy',
            # 'temperature': 1.0 * self('K'),
            'time': 'time',
            'x': 'length',
            'y': 'length',
            'z': 'length',
            'xyz_x': 'length',
            'xyz_y': 'length',
            'xyz_z': 'length',
            'position': 'length',
            'position_x': 'length',
            'position_y': 'length',
            'position_z': 'length',
            'dx': 'length',
            'mass': 'mass'
        }

    def __call__(self, *args, **kwargs):
        return self._ureg(*args, **kwargs).units

    def define(self, *args, **kwargs):
        self._ureg.define(*args, **kwargs)

    def set_base_units(self, ud, ul, ut):
        # self._base_units.update({
        self._base_units["density"] = ud * self("g / cm**3")
        self._base_units["velocity"] = (ul / ut) * self("cm / s")
        self._base_units["magnetic_field"] = sqrt(4.0 * pi * ud *
                                                  (ul / ut)**2) * self("G")
        self._base_units[
            "momentum"] = self._base_units["density"] * self._base_units["velocity"]
        self._base_units["acceleration"] = (ul / ut**2) * self("cm / s**2")
        self._base_units["energy"] = ud * ((ul / ut)**2) * self("erg / cm**3")
        self._base_units["time"] = ut * self("s")
        self._base_units["length"] = ul * self("cm")
        self._base_units[
            "mass"] = self._base_units["density"] * self._base_units["length"]**3
        self._base_units["temperature"] = self("K")

    def get(self, string):

        if string in self._ramses_units:
            return self._base_units[self._ramses_units[string]]
        return self._ureg('dimensionless')
        # return ramses_units.get(string, self('dimensionless')).units

        # if string in config.parameters['units']:
        #     return self(config.parameters['units'][string]).units
        # return ramses_units.get(string, self('dimensionless')).units


units = Units()
