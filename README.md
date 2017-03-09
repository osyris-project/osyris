# Plotting Ramses #

This is a small collection of python plotting scripts for RAMSES data. It is not meant to replace large projects such as Pymses, it's purpose is more to plot small 'quick and dirty' diagnostics while a simulation is running.

### Installation ###

You must first run 'f2py' on the fortran subroutine which reads in the RAMSES data:

```
#!bash

f2py -c read_ramses_data.f90 -m read_ramses_data
```

### From within ipython ###

```
#!python
import plotting_ramses as pp
mydata = pp.ramses_output(71)
pp.plot_histogram(mydata,'logrho','logT')
```

### From the terminal ###


```
#!bash

python make_figures.py 71

```


### Contributors ###

* Neil Vaytet (owner)