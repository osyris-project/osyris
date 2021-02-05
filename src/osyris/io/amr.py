from .loader import Loader
from .units import get_unit

class AmrLoader(Loader):

    def __init__(self, scale, units):

        self.initialized = True
        # AMR grid variables
	    length_unit = get_unit("x", units["ud"], units["ul"], units["ut"])
	    if scale is not None:
	        scale = units(scale)
	        scaling = (length_unit.to(scale) / scale).magnitude * scale
	    else:
	        scaling = length_unit

	    self.variables.update({"level": {"read": True, "type": "i", "buffer": None, "pieces": {}, "unit": 1.0*units.dimensionless},
	                     "cpu": {"read": True, "type": "i", "buffer": None, "pieces": {}, "unit": 1.0*units.dimensionless},
	                     "x": {"read": True, "type": "d", "buffer": None, "pieces": {}, "unit": scaling},
	                     "y": {"read": True, "type": "d", "buffer": None, "pieces": {}, "unit": scaling},
	                     "z": {"read": True, "type": "d", "buffer": None, "pieces": {}, "unit": scaling},
	                     "dx": {"read": True, "type": "d", "buffer": None, "pieces": {}, "unit": scaling}
	                     })