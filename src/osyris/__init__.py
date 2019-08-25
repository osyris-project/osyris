# Import the config from "/home/user/.osyris/config if it exists.
# If it doesn't, try to create one by copying the default from the source.
# If that fails, just load the default.
import os
import sys
conf_dir = os.path.join(os.path.expanduser("~"), ".osyris")
sys.path.append(conf_dir)
try:
    import config_osyris as config
except ImportError:
    from shutil import copyfile
    this_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(conf_dir):
        os.mkdir(conf_dir)
    copyfile(os.path.join(this_dir, "config.py"),
             os.path.join(conf_dir, "config_osyris.py"))
    try:
        import config_osyris as config
    except ImportError:
        from . import config

from .load_ramses import RamsesData
from .plot_histogram import plot_histogram
from .plot_slice import plot_slice
from .plot_column_density import plot_column_density
from .plot_3d import plot_volume, plot_quiver
from . import ism_physics
from .vtk import to_vtk
from .interpolate import interpolate
