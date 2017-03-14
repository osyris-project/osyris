from pylab import *
import plotting_ramses as pp

# Read arguments
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = 1

# Load data
mydata = pp.RamsesOutput(nout=nout,center="auto",scale="au")

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

# Create new fields
mydata.new_field(name="log_rho",values=log10(mydata.get_values("density")),unit="g/cm3",label="log(Density)")
mydata.new_field(name="log_T",values=log10(mydata.get_values("temperature")),unit="K",label="log(T)")
mydata.new_field(name="log_B",values=log10(mydata.get_values("B")),unit="G",label="log(B)")

# Density vs B field
pp.plot_histogram(mydata.get("log_rho"),mydata.get("log_B"),var_z=mydata.get("level"),axes=ax1,cmap="YlGnBu")
# Density vs Temperature
pp.plot_histogram(mydata.get("log_rho"),mydata.get("log_T"),var_z=mydata.get("level"),axes=ax2,cmap="YlGnBu")

#x,z density slice with B field
pp.plot_slice(mydata,"log_rho",direction="y",vec="B",dx=dx,dy=dy,axes=ax3,streamlines=True)
# x,y density slice with velocity
pp.plot_slice(mydata,"log_rho",direction="z",vec="velocity",dx=dx,dy=dy,axes=ax4)
# x,y temperature slice with velocity
pp.plot_slice(mydata,"log_T",direction="z",vec="velocity",dx=dx,dy=dy,axes=ax5,cmap='hot')
# z,z density slice with velocity
pp.plot_slice(mydata,"log_rho",direction="y",vec="velocity",dx=dx,dy=dy,axes=ax6)

fig.savefig("plots.pdf",bbox_inches="tight")
