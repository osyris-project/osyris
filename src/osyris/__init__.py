from .load_ramses import RamsesData
from .plot_histogram import plot_histogram
from .plot_slice import plot_slice
from .plot_column_density import plot_column_density
from .plot_3d import plot_volume
from . import ism_physics
from .vtk import to_vtk
from .interpolate import interpolate

# Import the config from "/home/user/.osyris/config if it exists.
# If it doesn't, try to create one by copying the default from the source.
# If that fails, just load the default.
from os.path import expanduser, dirname, abspath
this_dir = dirname(abspath(__file__))
from os import mkdir
conf_dir = expanduser("~") + "/.osyris"
import sys
sys.path.append(conf_dir)
try:
    import config
except ImportError:
    from shutil import copyfile
    mkdir(conf_dir)
    copyfile(this_dir + "/config.py", conf_dir + "/config.py")
    try:
        import config
    except ImportError:
        from . import config
    # try:
    #     from shutil import copyfile
    #     mkdir(conf_dir)
    #     copyfile("config.py", conf_dir + "/config.py")
    #     try:
    #         import config
    #     except ImportError:
    #         from . import config
    # except OSError:
    #     print("imported this")
    #     from . import config
