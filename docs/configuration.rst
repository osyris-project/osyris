.. _configuration:

Configuration
=============

Upon first import, ``osyris`` will create a configuration file located in
``/home/user/.osyris/config_osyris.py``.
This will then be loaded by ``osyris``, and allows easy configuration for the
default colormap, spatial scale, time unit etc.

Default values
--------------

- ``"nout"`` : 1 : The output number to load.
- ``"lmax"`` : 0 : The maximum level of cells that will be loaded
  (0 = no limit).
- ``"center"`` : None : The center around which the grid will be centered.
- ``"dx"`` : 0.0 : The ``x`` extent of the computational domain to be loaded.
- ``"dy"`` : 0.0 : The ``y`` extent of the computational domain to be loaded.
- ``"dz"`` : 0.0 : The ``z`` extent of the computational domain to be loaded.
- ``"scale"`` : "au" : The scale for the ``x, y, z`` coordinates.
- ``"time_unit"`` : "kyr" : The time unit (for labels, titles, etc.).
- ``"verbose"`` : False : Print data information when loaded if ``True``.
- ``"path"`` : "" : Path to the directory where to read outputs from,
  if different from the current directory.
- ``"variables"`` : [] : List of variables to be loaded.
- ``"colormap"`` : "viridis" : Color map for slices, images, etc.


Physical constants
------------------
- ``"cm"`` : 1.0 : Centimeter.
- ``"au"`` : 1.495980e+13 : Astronomical unit.
- ``"pc"`` : 3.085678e+18 : Parsec.
- ``"s"`` : 1.0 : Second.
- ``"yr"`` : 365.25*86400.0 : Year.
- ``"kyr"`` : 365.25*86400.0*1000.0 : Kilo-year.
- ``"msun"`` : 1.9889e33 : Solar mass.
- ``"a_r"`` : 7.56591469318689378e-015 : Radiation constant.
- ``"c"`` : 2.9979250e+10 : Speed of light.


Additional colormaps
--------------------

- ``"osyris"`` :
- ``"osyris2"`` :
- ``"osyris3"`` :
- ``"osyris4"`` :
- ``"osyris5"`` :

Additional variables
--------------------

Define variables, using either operations on existing variables or new values,
that will be created every time an output is loaded.
