.. _load-ramses-data:

Reading a RAMSES output
=======================

Navigate to the directory containing the data of your simulation and open an
``ipython`` console:

.. code-block:: sh

   $ cd /path/to/my/ramses/data
   ipython

Import the ``osyris`` library and load the output of your choice (this will be
output number 71 in this example).
**The data loader searches for a ``hydro_file_descriptor.txt`` inside the output
directory to get the variable names, so make sure your version of RAMSES
supports this.
If it doesn't, you can edit the ``var_names`` list in the `config.py`
configuration file, under ``default_values`` to match your data structure.
By default it will try to guess by itself which are the variables to read,
but this will almost certainly fail without editing it!**

(**Note:** you can download the sample data used in this tutorial
[here](http://www.nbi.dk/~nvaytet/osyris/ramses_sample_data.tar.gz).)

.. code-block:: python

   In [1]: import osyris
   In [2]: mydata = osyris.RamsesData(71,scale="au",verbose=True)
   ============================================
   Processing 60 files in output_00071
    10% : read     421904 cells
    20% : read     659832 cells
    30% : read     745768 cells
    40% : read    1208104 cells
    50% : read    1404744 cells
    60% : read    1826640 cells
    70% : read    2064568 cells
    80% : read    2150496 cells
    90% : read    2598376 cells
   Total number of cells loaded: 2809480
   Read 1 sink particles
   Generating data structure... please wait
   Memory used: 854.08 Mb
   output_00071 successfully loaded
   ---------------------------------------------------------------------------------------
   H0: 1.0
   aexp: 1.0
   boxlen: 0.0535934835282699
   boxsize: 1.650679292670713e+17
   boxsize_scaled: 11034.100005820352
   center: None
   dtnew: [0.0 ... 1.69759664423e-313]
   dtold: [0.0 ... 0.0]
   dx_load: 0.0
   dy_load: 0.0
   dz_load: 0.0
   eos: 1
   gamma: 1.666667
   infile: output_00071
   ir_cloud: 4
   leafs: [[   3152    3153    3154 ... 2809477 2809478 2809479] ... [   3152    3153    3154 ... 2809477 2809478 2809479]]
   levelmax: 29
   levelmax_active: 17.0
   levelmin: 6
   lmax: 0
   mu_gas: 2.31
   ncells: 2809480
   ncpu: 60
   ndim: 3
   ngridmax: 400000
   ngrp: 1
   nout: 71
   npart_tot: 0
   nsinks: 1
   nstep_coarse: 700
   nvar_hydro: 17
   omega_b: 0.0
   omega_k: 0.0
   omega_l: 0.0
   omega_m: 1.0
   ordering type: hilbert
   path:
   scale: au
   time: 899016308300.144
   unit_d: 3.8346e-24
   unit_l: 3.08e+18
   unit_t: 1977320409468800.0
   variables: []
   xc: 5517.050002910176
   yc: 5517.050002910176
   zc: 5517.050002910176
   ---------------------------------------------------------------------------------------
   The variables are:
   Name               Type   Group Unit      Min                    Max
   B                  vector hydro [G      ] 7.0556866300589425e-06 18.122346368805925
   B_left             vector hydro [G      ] 7.049973426689696e-06  18.19666750359762
   B_right            vector hydro [G      ] 7.049973426372603e-06  18.196664964960195
   cpu                scalar amr   [       ] 1.0                    60.0
   density            scalar hydro [g/cm3  ] 1.5375905866344673e-20 2.6285126781491883e-09
   dx                 scalar amr   [au     ] 0.08418350224167138    172.407812590943
   grav_acceleration  vector grav  [cm/s2  ] 2.803138592806718e-10  0.0050114715975282265
   grav_potential     scalar grav  [       ] -1523069.2164004715    4188.759353776809
   leaf               scalar amr   [       ] 1.0                    1.0
   level              scalar amr   [       ] 6.0                    17.0
   log_B              scalar hydro [G      ] -5.151460716121474     1.258214426732307
   log_T              scalar hydro [K      ] 0.9779806731254507     2.8482524904907063
   log_m              scalar hydro [Msun   ] -9.509941029873582     -4.548100746783424
   log_r              scalar amr   [au     ] -1.1372423788049735    3.9734281229644557
   log_rho            scalar hydro [g/cm3  ] -19.813159288484858    -8.580289923900661
   mass               scalar hydro [Msun   ] 3.090715072856166e-10  2.830735251207072e-05
   passive_scalar_1   scalar hydro [       ] 0.0                    0.0
   passive_scalar_2   scalar hydro [       ] 0.0                    0.0
   passive_scalar_3   scalar hydro [       ] 0.0                    0.0
   passive_scalar_4   scalar hydro [       ] 209.4555017278068      24103.182593448306
   r                  scalar amr   [au     ] 0.07290505151880808    9406.501367423783
   radiative_energy_1 scalar hydro [erg/cm3] 6.247691684509019e-11  0.0018699559894006896
   temperature        scalar hydro [K      ] 9.505624909986         705.1028831330501
   thermal_pressure   scalar hydro [erg/cm3] 5.236332931938116e-12  102.4803877152337
   velocity           vector hydro [cm/s   ] 157.56222171237246     320341.487894837
   x                  scalar amr   [au     ] -5430.846096614704     5430.846096614703
   y                  scalar amr   [au     ] -5430.846096614704     5430.846096614703
   z                  scalar amr   [au     ] -5430.846096614704     5430.846096614703
   ============================================

In the call to ``RamsesData``, the first argument is the output number.
**Note:** you can use ``-1`` to select the last output in the directory.
The second argument is the spatial scale you want to convert distances to.
Possible choices are ``"cm"``, ``"au"`` or ``"pc"``.
If you add ``verbose=True`` to the argument list, it will also print out some
information about the data (the variables names, their minimum and maximum
values, etc.). ``osyris`` tries to guess the units of each variable field
according to its name. This is done by the ``get_units()`` function and can
easily be modified if you have non-standard variables.