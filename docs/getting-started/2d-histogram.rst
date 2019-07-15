.. _2d-histogram:

2D histogram
============

We now wish to plot a 2d histogram of the logarithm of density `"log_rho"` versus logarithm of gas temperature `"log_T"` for all the cells inside the computational domain. We also use a logarithmic colormap which represents the cell density in the (rho,T) plane.
```
#!python

In [3]: osiris.plot_histogram(mydata.log_rho,mydata.log_T,scalar_args={"cmap":"log"})

```
This creates a figure which looks like

![figure_1-25.png](https://bitbucket.org/repo/jq5boX/images/292924844-figure_1-25.png)

You can also save the figure to file directly by specifying the argument `fname="rho_T_histogram.pdf"` in the call.
