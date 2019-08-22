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

.. raw:: html

   <ul>
   <li>"osyris" :
   <svg width="414" height="24" >
     <linearGradient id="osyris" x1="0" x2="1" y1="0" y2="0">
       <stop style="stop-color:#2b3c4e;stop-opacity:1" offset="0%"/>
       <stop style="stop-color:#249593;stop-opacity:1" offset="33.3%"/>
       <stop style="stop-color:#db6a6c;stop-opacity:1" offset="66.7%"/>
       <stop style="stop-color:#ffffff;stop-opacity:1" offset="100%"/>
     </linearGradient>
     <rect x="10" y="2" width="400" height="20" style="fill:url(#osyris);stroke-width:2;stroke:#000000"/>
   </svg>
   <li>"osyris2" :
   <svg width="404" height="24" >
     <linearGradient id="osyris2" x1="0" x2="1" y1="0" y2="0">
       <stop style="stop-color:#2b3c4e;stop-opacity:1" offset="0%"/>
       <stop style="stop-color:#249593;stop-opacity:1" offset="25%"/>
       <stop style="stop-color:#ffffff;stop-opacity:1" offset="50%"/>
       <stop style="stop-color:#db6a6c;stop-opacity:1" offset="75%"/>
       <stop style="stop-color:#9e4d4e;stop-opacity:1" offset="100%"/>
     </linearGradient>
     <rect x="2" y="2" width="400" height="20" style="fill:url(#osyris2);stroke-width:2;stroke:#000000"/>
   </svg>
   <li>"osyris3" :
   <svg width="404" height="24" >
     <linearGradient id="osyris3" x1="0" x2="1" y1="0" y2="0">
       <stop style="stop-color:#3d3d6b;stop-opacity:1" offset="0%"/>
       <stop style="stop-color:#2a7b9b;stop-opacity:1" offset="10%"/>
       <stop style="stop-color:#00baad;stop-opacity:1" offset="20%"/>
       <stop style="stop-color:#57c785;stop-opacity:1" offset="30%"/>
       <stop style="stop-color:#add45c;stop-opacity:1" offset="40%"/>
       <stop style="stop-color:#ffc300;stop-opacity:1" offset="50%"/>
       <stop style="stop-color:#ff8d1a;stop-opacity:1" offset="60%"/>
       <stop style="stop-color:#ff5733;stop-opacity:1" offset="70%"/>
       <stop style="stop-color:#c70039;stop-opacity:1" offset="80%"/>
       <stop style="stop-color:#900c3f;stop-opacity:1" offset="90%"/>
       <stop style="stop-color:#511849;stop-opacity:1" offset="100%"/>
     </linearGradient>
     <rect x="2" y="2" width="400" height="20" style="fill:url(#osyris3);stroke-width:2;stroke:#000000"/>
   </svg>
   <li>"osyris4" :
   <svg width="404" height="24" >
     <linearGradient id="osyris4" x1="0" x2="1" y1="0" y2="0">
       <stop style="stop-color:#000000;stop-opacity:1" offset="0%"/>
       <stop style="stop-color:#ff5b00;stop-opacity:1" offset="14.3%"/>
       <stop style="stop-color:#ffff00;stop-opacity:1" offset="28.6%"/>
       <stop style="stop-color:#00ff00;stop-opacity:1" offset="42.9%"/>
       <stop style="stop-color:#2bc184;stop-opacity:1" offset="57.1%"/>
       <stop style="stop-color:#3d3d6b;stop-opacity:1" offset="71.4%"/>
       <stop style="stop-color:#ffffff;stop-opacity:1" offset="85.7%"/>
       <stop style="stop-color:#0000ff;stop-opacity:1" offset="100%"/>
     </linearGradient>
     <rect x="2" y="2" width="400" height="20" style="fill:url(#osyris4);stroke-width:2;stroke:#000000"/>
   </svg>
   <li>"osyris5" :
   <svg width="404" height="24" >
     <linearGradient id="osyris5" x1="0" x2="1" y1="0" y2="0">
       <stop style="stop-color:#0000ff;stop-opacity:1" offset="0%"/>
       <stop style="stop-color:#00ffff;stop-opacity:1" offset="16.7%"/>
       <stop style="stop-color:#00ff00;stop-opacity:1" offset="33.3%"/>
       <stop style="stop-color:#ffff00;stop-opacity:1" offset="50%"/>
       <stop style="stop-color:#ff0000;stop-opacity:1" offset="66.7%"/>
       <stop style="stop-color:#000000;stop-opacity:1" offset="83.3%"/>
       <stop style="stop-color:#ffffff;stop-opacity:1" offset="100%"/>
     </linearGradient>
     <rect x="2" y="2" width="400" height="20" style="fill:url(#osyris5);stroke-width:2;stroke:#000000"/>
   </svg>
   </ul>

Additional variables
--------------------

Define variables, using either operations on existing variables or new values,
that will be created every time an output is loaded.

For instance, to create the logarithm of gas density, we use

.. code-block:: python

   holder.new_field(name="log_rho", operation="np.log10(density)",
                    unit="g/cm3", label="log(Density)", verbose=False)
