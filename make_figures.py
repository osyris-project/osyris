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
ratio = 1.4
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)
ax1  = fig.add_subplot(211)
ax2  = fig.add_subplot(212)

# Density - B field
#pp.plot_histogram(mydata,"logrho","logB",fname="brho.pdf")
pp.plot_histogram(mydata,"logrho","logB",axes=ax1)

# Density - Temperature
#pp.plot_histogram(mydata,"logrho","logT",fname="trho.pdf")
pp.plot_histogram(mydata,"logrho","logT",axes=ax2)

# x,y density slice
#plot_slice(data,"x","y","rho",fname="rhoxy.pdf",zlog=True,xmin=-50.0,xmax=50.0,ymin=-50.0,ymax=50.0)

fig.savefig("plots.pdf",bbox_inches="tight")
