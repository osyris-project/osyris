![demo.png](https://bitbucket.org/repo/jq5boX/images/2936418214-demo.png)

# Osiris #

A python visualization utility for RAMSES data. It's purpose is to plot quick diagnostics while a simulation is running, and also produce publication grade figures.

### Installation ###

You will need matplotlib installed on your system. Clone the Osiris repository and append its location to your PYTHONPATH.

### From within ipython ###

```
#!python
import osiris as pp
mydata = pp.RamsesData(71,scale="au")
mydata.plot_slice("log_rho",direction="z",vec="velocity",dx=100)
```

### Demo ###

You can download the sample data [here](http://www.nbi.dk/~nvaytet/osiris/ramses_sample_data.tar.gz).

```
#!python

import matplotlib.pyplot as plt
import numpy as np
import osiris as pp

# Change default time unit to kyr
pp.conf.default_values["time_unit"]="kyr"

# Load data
mydata = pp.RamsesData(nout=71,center="max:density",scale="au")

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

# Density vs B field with AMR level contours
mydata.plot_histogram("log_rho","log_B",var_c="level",axes=ax1,cmap="YlGnBu",contour_args={"fmt":"%i"})

# Create new field with log of velocity
mydata.new_field(name="log_vel",operation="np.log10(np.sqrt(velocity_x**2+velocity_y**2+velocity_z**2))",unit="cm/s",label="log(Velocity)")

# Density vs log_vel
mydata.plot_histogram("log_rho","log_vel","log_T",axes=ax2,cmap="gnuplot",scatter=True,outline=True,scatter_args={"iskip":100})

#x,z density slice with B field streamlines
mydata.plot_slice("log_rho",direction="y",stream="B",dx=100,axes=ax3)
# x,y density slice with velocity vectors in color
mydata.plot_slice("log_rho",direction="z",vec="velocity",dx=100,axes=ax4,vec_args={"cmap":"seismic","vskip":4})
# x,y temperature slice with velocity vectors
mydata.plot_slice("log_T",direction="z",vec="velocity",dx=100,axes=ax5,cmap="hot")

# Now update values with later snapshot
mydata.update_values(201)
# Re-plot x,y density slice with velocity vectors in color
mydata.plot_slice("log_rho",direction="auto:top",vec="velocity",dx=100,axes=ax6)

fig.savefig("demo.pdf",bbox_inches="tight")

```


### Contributors ###

* Neil Vaytet (StarPlan/NBI)
* Tommaso Grassi (StarPlan/NBI)
* Matthias Gonz�lez (CEA Saclay)