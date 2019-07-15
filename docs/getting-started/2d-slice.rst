.. _2d-slice:

2D slice
========

To take a slice through the data, simply type

.. code-block:: python

   osiris.plot_slice(mydata.log_rho, direction="z", vec=mydata.velocity, dx=100)

where the first argument is the variable to display, ``direction`` is the normal
to the slice plane, ``vec`` is the (optional) variable to be used to plot
vectors, and ``dx`` is the extent of the slice.

**Note:** the units of ``dx`` are consistent with the ``scale`` specified when
reading in the snapshot using ``RamsesOutput``.

This should produce the following figure

![figure_1-26.png](https://bitbucket.org/repo/jq5boX/images/2245240735-figure_1-26.png)
