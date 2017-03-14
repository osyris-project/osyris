![demo.png](https://bitbucket.org/repo/jq5boX/images/1336351696-demo.png)
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

You can also plot a small collection of plots to a pdf file:
```
#!bash

python make_figures.py 71

```

### Short example ###


```
#!python

from pylab import *
import plotting_ramses as pp

# Load data: select only cells 100 au on each side of a center located at (0.5,0.5,0.5)
mydata = pp.RamsesOutput(nout=71,center=[0.5,0.5,0.5],scale="au",dx=200,dy=200,dz=200)

# Create figure
fig = matplotlib.pyplot.figure()
ratio = 0.5
sizex = 20.0
fig.set_size_inches(sizex,ratio*sizex)
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

# Create new fields
mydata.new_field(name="log_rho",values=log10(mydata.get_values("rho")),unit="g/cm3",label="log(Density)")
mydata.new_field(name="log_T",values=log10(mydata.get_values("T")),unit="K",label="log(T)")
mydata.new_field(name="log_B",values=log10(mydata.get_values("B")),unit="G",label="log(B)")

# Histogram Density vs B field, with overlayed AMR level contours
pp.plot_histogram(mydata.get("log_rho"),mydata.get("log_B"),var_z=mydata.get("level"),axes=ax1,cmap="YlGnBu")

# You can also feed simple data arrays to the plot_histogram routine, it does not have to be the "mydata" fields
# Note that in this case you will not have labels on the plot axes
# Histogram Density vs Temperature
pp.plot_histogram(log10(mydata.get_values("rho")),log10(mydata.get_values("T")),var_z=mydata.get_values("level"),axes=ax2,cmap="YlGnBu")

# Define size of the slices
dx = 100
dy = 100

#x,z density slice with B field streamlines
pp.plot_slice(mydata,"log_rho",direction="y",vec="B",dx=dx,dy=dy,axes=ax3,streamlines=True)
# x,y density slice with velocity vectors
pp.plot_slice(mydata,"log_rho",direction="z",vec="vel",dx=dx,dy=dy,axes=ax4)
# x,y temperature slice with velocity vectors
pp.plot_slice(mydata,"log_T",direction="z",vec="vel",dx=dx,dy=dy,axes=ax5,cmap='hot')
# x,z density slice with velocity vectors
pp.plot_slice(mydata,"log_rho",direction="y",vec="vel",dx=dx,dy=dy,axes=ax6)

fig.savefig("plots.pdf",bbox_inches="tight")
```


### Contributors ###

* Neil Vaytet (owner)