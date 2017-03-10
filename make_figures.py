from pylab import *
import plotting_ramses as pp

# Read arguments
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = 1

# Load data
mydata=pp.ramses_output(nout=nout)

# Create figure
fig = matplotlib.pyplot.figure()
ratio = 0.5
sizex = 20.0
fig.set_size_inches(sizex,ratio*sizex)
ax1  = fig.add_subplot(231)
ax2  = fig.add_subplot(232)
ax3  = fig.add_subplot(234)
ax4  = fig.add_subplot(235)
ax5  = fig.add_subplot(236)

# Density - B field
pp.plot_histogram(mydata,"logrho","logB",axes=ax1,cmap="YlGnBu")

# Density - Temperature
pp.plot_histogram(mydata,"logrho","logT",axes=ax2,cmap="YlGnBu")

# x,y density slice
pp.plot_slice(mydata,'z_au','logrho',dx=100,dy=100,axes=ax3)

# x,y density slice
pp.plot_slice(mydata,'y_au','logrho',dx=100,dy=100,axes=ax4)

# x,y density slice
pp.plot_slice(mydata,'x_au','logrho',dx=100,dy=100,axes=ax5)


fig.savefig("plots.pdf",bbox_inches="tight")
