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

# Physical constants
au = 1.495980e+13

# Load data: select only cells 100 au on each side of a center located at (0.5,0.5,0.5)
mydata = pp.RamsesOutput(nout=71,center=[0.5,0.5,0.5],scale=au,dx=200,dy=200,dz=200)

# Create figure
fig = matplotlib.pyplot.figure()
ratio = 0.5
sizex = 20.0
fig.set_size_inches(sizex,ratio*sizex)
ax1  = fig.add_subplot(231)
ax2  = fig.add_subplot(232)
ax3  = fig.add_subplot(233)
ax4  = fig.add_subplot(234)
ax5  = fig.add_subplot(235)
ax6  = fig.add_subplot(236)

dx = 100
dy = 100

# Histogram Density vs B field
pp.plot_histogram(log10(mydata.rho),log10(mydata.B[:,3]),axes=ax1,cmap="YlGnBu")
# Histogram Density vs Temperature
pp.plot_histogram(log10(mydata.rho),log10(mydata.T),axes=ax2,cmap="YlGnBu")
# x,z density slice with B field streamlines
pp.plot_slice(mydata.x,log10(mydata.rho),direction=1,vec=mydata.B,dx=dx,dy=dy,axes=ax3,streamlines=True)
# x,y density slice with velocity vectors
pp.plot_slice(mydata.x,log10(mydata.rho),direction=2,vec=mydata.vel,dx=dx,dy=dy,axes=ax4)
# y,z temperature slice with velocity vectors
pp.plot_slice(mydata.x,log10(mydata.T),direction=0,vec=mydata.vel,dx=dx,dy=dy,axes=ax5,cmap='hot')
# x,z density slice with velocity vectors
pp.plot_slice(mydata.x,log10(mydata.rho),direction=1,vec=mydata.vel,dx=dx,dy=dy,axes=ax6)

fig.savefig("plots.pdf",bbox_inches="tight")
```


### Contributors ###

* Neil Vaytet (owner)