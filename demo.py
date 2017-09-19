#=======================================================================================
#This file is part of OSIRIS.

#OSIRIS is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#OSIRIS is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with OSIRIS.  If not, see <http://www.gnu.org/licenses/>.
#=======================================================================================

import matplotlib.pyplot as plt
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
mydata.plot_histogram("log_rho","log_B",axes=ax1,cmap="log,YlGnBu")
mydata.plot_histogram("log_rho","log_B",var_z="level",contour=True,axes=ax1,contour_args={"fmt":"%i","label":True,"colors":"k","cmap":None,"levels":range(5,20)},cbar=False,zmin=6,zmax=16)

# Create new field with log of velocity
mydata.new_field(name="log_vel",operation="np.log10(np.sqrt(velocity_x**2+velocity_y**2+velocity_z**2))",unit="cm/s",label="log(Velocity)")

# Density vs log_vel
mydata.plot_histogram("log_rho","log_vel","log_T",axes=ax2,cmap="gnuplot",scatter=True,outline=True,scatter_args={"iskip":100})

#x,z density slice with B field streamlines
mydata.plot_slice("density",direction="y",stream="B",dx=100,axes=ax3,cmap="log")
# x,y density slice with velocity vectors in color
mydata.plot_slice("log_rho",direction="z",vec="velocity",dx=100,axes=ax4,vec_args={"cmap":"seismic","vskip":4})
# x,y temperature slice with velocity vectors
mydata.plot_slice("log_T",direction="z",vec="velocity",dx=100,axes=ax5,cmap="hot")
mydata.plot_slice("level",direction="z",dx=100,axes=ax5,contour=True,contour_args={"fmt":"%i","label":True,"colors":"w","cmap":None,"levels":range(9,17)},cbar=False)

# Now update values with later snapshot
mydata.update_values(201)
# Re-plot x,y density slice with velocity vectors in color
mydata.plot_slice("log_rho",direction="auto:top",vec="velocity",dx=100,axes=ax6)

fig.savefig("demo.pdf",bbox_inches="tight")
