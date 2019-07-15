.. _load-ramses-data:

Reading a RAMSES output
=======================

Navigate to the directory containing the data of your simulation and open an `ipython` console:
```
$ cd /path/to/my/ramses/data
ipython
```
Import the `osiris` library and load the output of your choice (this will be output number 71 in this example). **The data loader searches for a `hydro_file_descriptor.txt` inside the output directory to get the variable names, so make sure your version of RAMSES supports this. If it doesn't, you can edit the `"var_names` list in the `config_osiris.py` configuration file, under `default_values` to match your data structure. By default it will try to guess by itself which are the variables to read, but this will almost certainly fail without editing it!**

(**Note:** you can download the sample data used in this tutorial [here](http://www.nbi.dk/~nvaytet/osiris/ramses_sample_data.tar.gz).)

.. code-block:: python

   In [1]: import osiris
   In [2]: mydata = osiris.RamsesData(71,scale="au",verbose=True)
   ============================================
   Processing 60 files in output_00071
    10% : read     369140 cells
    20% : read     577355 cells
    30% : read     652555 cells
    40% : read    1057064 cells
    50% : read    1229148 cells
    60% : read    1598288 cells
    70% : read    1806511 cells
    80% : read    1881703 cells
    90% : read    2273577 cells
   Total number of cells loaded: 2458296
   Generating data structure... please wait
   Memory used: 884.99 Mb
   output_00071 successfully loaded
   --------------------------------------------
   H0: 1.0
   aexp: 1.0
   boxlen: 0.0535934835283
   boxsize: 1.65067929267e+17
   boxsize_scaled: 11034.1000058
   center: max:density
   dtnew: [0.0 ... 0.0]
   dtold: [0.0 ... 0.0]
   dx_load: 0.0
   dy_load: 0.0
   dz_load: 0.0
   eos: 1
   gamma: 1.666667
   infile: output_00071
   ir_cloud: 4
   levelmax: 29
   levelmin: 6
   lmax: 0
   mu_gas: 2.31
   ncells: 2458296
   ncpu: 60
   ndim: 3
   ngridmax: 400000
   ngrp: 1
   nout: 71
   nsinks: 0
   nstep_coarse: 700
   nvar_hydro: 17
   omega_b: 0.0
   omega_k: 0.0
   omega_l: 0.0
   omega_m: 1.0
   ordering type: hilbert
   path:
   scale: au
   time: 8.990163083e+11
   unit_d: 3.8346e-24
   unit_l: 3.08e+18
   unit_t: 1.97732040947e+15
   variables: []
   xc: 5517.09209466
   yc: 5517.00791116
   zc: 5517.00791116
   --------------------------------------------
   The variables are:
   Name               Type   Group  Unit      Min               Max
   B                  vector hydro [G      ] 7.05568663006e-06 18.1223463688
   B_left             vector hydro [G      ] 7.04997342669e-06 18.1966675036
   B_right            vector hydro [G      ] 7.04997342637e-06 18.196664965
   cpu                scalar amr   [       ] 1.0               60.0
   density            scalar hydro [g/cm3  ] 1.53759058663e-20 2.62851267815e-09
   dx                 scalar amr   [au     ] 0.0841835022417   172.407812591
   dx_box             scalar amr   [       ] 7.62939453125e-06 0.015625
   dx_raw             scalar amr   [au     ] 1.25936835683e+12 2.5791863948e+15
   grav_acceleration  vector grav  [       ] 355.834417609     6361633641.12
   grav_potential     scalar grav  [       ] -1523069.2164     4188.75935378
   level              scalar amr   [       ] 6.0               17.0
   log_B              scalar hydro [G      ] -5.15146071612    1.25821442673
   log_T              scalar hydro [K      ] 0.977980673125    2.84825249049
   log_m              scalar hydro [Msun   ] -9.50994102987    -4.54810074678
   log_r              scalar amr   [au     ] -inf              3.97343148895
   log_rho            scalar hydro [g/cm3  ] -19.8131592885    -8.5802899239
   mass               scalar hydro [Msun   ] 3.09071507286e-10 2.83073525121e-05
   passive_scalar_1   scalar hydro [       ] 0.0               0.0
   passive_scalar_2   scalar hydro [       ] 0.0               0.0
   passive_scalar_3   scalar hydro [       ] 0.0               0.0
   passive_scalar_4   scalar hydro [       ] 209.455501728     24103.1825934
   r                  scalar amr   [au     ] 0.0               9406.57427248
   radiative_energy_1 scalar hydro [erg/cm3] 6.24769168451e-11 0.0018699559894
   temperature        scalar hydro [K      ] 9.50562490999     705.102883133
   thermal_pressure   scalar hydro [erg/cm3] 5.23633293194e-12 102.480387715
   velocity           vector hydro [cm/s   ] 157.562221712     320341.487895
   x                  scalar amr   [au     ] -5430.88818837    5430.80400486
   x_box              scalar amr   [       ] 0.0078125         0.9921875
   x_raw              scalar amr   [au     ] 1.2895931974e+15  1.6377833607e+17
   y                  scalar amr   [au     ] -5430.80400486    5430.88818837
   y_box              scalar amr   [       ] 0.0078125         0.9921875
   y_raw              scalar amr   [au     ] 1.2895931974e+15  1.6377833607e+17
   z                  scalar amr   [au     ] -5430.80400486    5430.88818837
   z_box              scalar amr   [       ] 0.0078125         0.9921875
   z_raw              scalar amr   [au     ] 1.2895931974e+15  1.6377833607e+17
   ============================================

In the call to `RamsesData`, the first argument is the output number. **Note:** you can use `-1` to select the last output in the directory. The second argument is the spatial scale you want to convert distances to. Possible choices are `"cm"`, `"au"` or `"pc"`.
If you add `verbose=True` to the argument list, it will also print out some information about the data (the variables names, their minimum and maximum values, etc.). `osiris` tries to guess the units of each variable field according to its name. This is done by the `get_units()` function and can easily be modified if you have non-standard variables.