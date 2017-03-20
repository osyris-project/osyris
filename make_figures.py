import matplotlib.pyplot as plt
import numpy as np
import plotting_ramses as pp
import sys

# Read arguments
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = 1

# Load data
mydata = pp.RamsesOutput(nout=nout,center="auto",scale="au")

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

dx = 100

# Density vs B field
mydata.plot_histogram("log_rho","log_B",var_z="level",axes=ax1,cmap="YlGnBu")
# Density vs Temperature
mydata.plot_histogram("log_rho","log_T",axes=ax2,cmap="YlGnBu")

#x,z density slice with B field
mydata.plot_slice("log_rho",direction="y",vec="B",dx=dx,axes=ax3,streamlines=True)
# x,y density slice with velocity
mydata.plot_slice("log_rho",direction="z",vec="velocity",dx=dx,axes=ax4)
# x,y temperature slice with velocity
mydata.plot_slice("log_T",direction="z",vec="velocity",dx=dx,axes=ax5,cmap="hot")
# x,z density slice with velocity
mydata.plot_slice("log_rho",direction="y",vec="velocity",dx=dx,axes=ax6)

fig.savefig("plots.pdf",bbox_inches="tight")
