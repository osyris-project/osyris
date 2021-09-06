import numpy as np
from .. import units


def get_unit(string, ud, ul, ut):
    ramses_units = {
        "density": ud * (units.g / (units.cm**3)),
        "velocity": (ul / ut) * (units.cm / units.s),
        "part_velocity": (ul / ut) * (units.cm / units.s),
        "momentum": (ud * ul / ut) * (units.g / (units.cm**2) / units.s),
        "B_": np.sqrt(4.0 * np.pi * ud * (ul / ut)**2) * units.G,
        "part_tracer_b": np.sqrt(4.0 * np.pi * ud * (ul / ut)**2) * units.G,
        "acceleration": (ul / ut**2) * (units.cm / (units.s**2)),
        "thermal_pressure": ud * ((ul / ut)**2) * (units.erg / (units.cm**3)),
        "energy": ud * ((ul / ut)**2) * (units.erg / (units.cm**3)),
        "temperature": 1.0 * units.K,
        "time": ut * units.s
    }
    ramses_units.update(dict.fromkeys(['x', 'y', 'z', 'dx'], (ul * units.cm)))

    if string in ramses_units:
        return ramses_units[string]

    for key in ramses_units:
        if string.startswith(key):
            return ramses_units[key]

    for key in ramses_units:
        if key in string:
            return ramses_units[key]

    return 1.0 * units.dimensionless
