# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)
"""
Define default values so that you don't have to specify them every time.
"""
parameters = {
    "scale": "au",
    "path": None,
    "select": None,
    "cmap": "viridis",
    "render_mode": "pcolormesh"
}


def additional_units(ureg):
    """
    Define additional useful units and constants
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
        data["B_field"] = 0.5 * (data["B_left"] + data["B_right"])
    except KeyError:
        pass

    # Mass
    try:
        data["mass"] = data["density"] * data["dx"]**3
        data["mass"].to("msun")
    except KeyError:
        pass
