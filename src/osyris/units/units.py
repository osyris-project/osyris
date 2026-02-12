# SPDX-License-Identifier: BSD-3-Clause

from pint import Quantity, Unit, UnitRegistry, compat


class Units:
    def __init__(self):
        self._ureg = UnitRegistry(system="cgs")

        # Define some useful constants
        self._ureg.define("bolometric_luminosity = 3.0128e+28 * W = L_bol0")
        self._ureg.define("solar_luminosity = 3.828e+26 * W = L_sun = L_sol")
        self._ureg.define("earth_mass = 5.97216787e+27 * g = M_earth")
        self._ureg.define("jupiter_mass = 1.8981246e+30 * g = M_jup")
        self._ureg.define("solar_mass = 1.9889e+33 * g = M_sun = M_sol")
        self._ureg.define("earth_radius = 6.3781e+08 * cm = R_earth")
        self._ureg.define("jupiter_radius = 7.1492e+09 * cm = R_jup")
        self._ureg.define("solar_radius = 6.957e+10 * cm = R_sun = R_sol")
        self._ureg.define(
            "radiation_constant = 7.56591469318689378e-015 * erg / cm^3 / K^4 = ar"
        )

    def __call__(self, arg):
        if isinstance(arg, Quantity):
            raise TypeError("Cannot create unit from a Quantity. Use `.units` instead.")
        if isinstance(arg, Unit):
            return arg
        return self._ureg(arg).units

    def define(self, *args, **kwargs):
        self._ureg.define(*args, **kwargs)


units = Units()

# Add compatibility for upcasting
compat.upcast_type_map["osyris.core.array.Array"] = None
compat.upcast_type_map["osyris.core.vector.Vector"] = None
