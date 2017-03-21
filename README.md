![demo.png](https://bitbucket.org/repo/jq5boX/images/1336351696-demo.png)
# Osiris #

This is a small collection of python plotting scripts for RAMSES data. It is not meant to replace large projects such as Pymses, it's purpose is more to plot small 'quick and dirty' diagnostics while a simulation is running.

### Installation ###

You will need matplotlib and f2py installed on your system.
Before plotting, you must first run 'f2py' on the fortran subroutine which reads in the RAMSES data:

```
#!bash

f2py -c read_ramses_data.f90 -m read_ramses_data
```

### From within ipython ###

```
#!python
import osiris as pp
mydata = pp.RamsesOutput(71)
mydata.plot_histogram("log_rho","log_B")
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

import matplotlib.pyplot as plt
import numpy as np
import osiris as pp

# Load data: select only cells 100 au on each side of a center located at (0.5,0.5,0.5)
mydata = pp.RamsesOutput(nout=71,center=[0.5,0.5,0.5],scale="au",dx=200,dy=200,dz=200)

# Create figure
fig = plt.figure()
ratio = 0.5
sizex = 20.0
fig.set_size_inches(sizex,ratio*sizex)
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

# Histogram Density vs B field, with overlayed AMR level contours
mydata.plot_histogram("log_rho","log_B",var_z="level",axes=ax1,cmap="YlGnBu")

# Create new field
mydata.new_field(name="log_vel",operation="np.log10(np.sqrt(velocity_x**2+velocity_y**2+velocity_z**2))",unit="cm/s",label="log(Velocity)")

# Histogram Density vs Velocity
mydata.plot_histogram("log_rho","log_vel",axes=ax2,cmap="YlGnBu")

# Define size of the slices
dx = 100
dy = 100

#x,z density slice with B field streamlines
mydata.plot_slice("log_rho",direction="y",vec="B",dx=dx,dy=dy,axes=ax3,streamlines=True)
# x,y density slice with velocity vectors
mydata.plot_slice("log_rho",direction="z",vec="velocity",dx=dx,dy=dy,axes=ax4)
# x,y temperature slice with velocity vectors
mydata.plot_slice("log_T",direction="z",vec="velocity",dx=dx,dy=dy,axes=ax5,cmap='hot')
# x,z density slice with velocity vectors
mydata.plot_slice("log_rho",direction="y",vec="velocity",dx=dx,dy=dy,axes=ax6)

fig.savefig("plots.pdf",bbox_inches="tight")
```


### Contributors ###

* Neil Vaytet (StarPlan/NBI)
* Tommaso Grassi (StarPlan/NBI)