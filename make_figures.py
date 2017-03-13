from pylab import *
import plotting_ramses as pp

# Read arguments
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = 1

# Physical constants
au = 1.495980e+13

# Load data
mydata = pp.RamsesOutput(nout=nout,center="auto",scale=au)

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

# Density - B field
pp.plot_histogram(log10(mydata.rho),log10(mydata.B[:,3]),axes=ax1,cmap="YlGnBu")
# Density - Temperature
pp.plot_histogram(log10(mydata.rho),log10(mydata.T),axes=ax2,cmap="YlGnBu")
# x,y density slice
pp.plot_slice(mydata.x,log10(mydata.rho),direction=1,vec=mydata.B,dx=dx,dy=dy,axes=ax3,streamlines=True)
# x,y density slice
pp.plot_slice(mydata.x,log10(mydata.rho),direction=2,vec=mydata.vel,dx=dx,dy=dy,axes=ax4)
# x,z density slice
pp.plot_slice(mydata.x,log10(mydata.T),direction=0,vec=mydata.vel,dx=dx,dy=dy,axes=ax5,cmap='hot')
# y,z density slice
pp.plot_slice(mydata.x,log10(mydata.rho),direction=1,vec=mydata.vel,dx=dx,dy=dy,axes=ax6)

fig.savefig("plots.pdf",bbox_inches="tight")
