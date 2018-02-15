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

#=======================================================================================
# This file is a short introduction on how to use the OSIRIS visualization package
# For more demos, please see the wiki page at:
#     https://bitbucket.org/nvaytet/osiris/wiki/Demos
#=======================================================================================

import matplotlib.pyplot as plt
import osiris

# Change default time unit to kyr
osiris.conf.default_values["time_unit"]="kyr"

# Load data
mydata = osiris.RamsesData(nout=71,center="max:density",scale="au",verbose=True)

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
osiris.plot_histogram(mydata.log_rho,mydata.log_B,axes=ax1,scalar=True,scalar_args={"cmap":"log,YlGnBu"},contour=mydata.level,contour_args={"fmt":"%i","label":True,"colors":"k","cmap":None,"levels":range(5,20),"cbar":False})

# Create new field with log of velocity
mydata.new_field(name="log_vel",operation="np.log10(np.sqrt(velocity_x**2+velocity_y**2+velocity_z**2))",unit="cm/s",label="log(Velocity)")

# Density vs log_vel in scatter mode with a grey outline
osiris.plot_histogram(mydata.log_rho,mydata.log_vel,axes=ax2,scatter=mydata.log_T,scatter_args={"iskip":100,"cmap":"gnuplot"},outline=True)

#x,z density slice with B field streamlines
osiris.plot_slice(mydata.density,direction="y",stream=mydata.B,dx=100,axes=ax3,scalar_args={"cmap":"log"})
# x,y density slice with velocity vectors in color
osiris.plot_slice(scalar=mydata.log_rho,direction="z",vec=mydata.velocity,dx=100,axes=ax4,vec_args={"cmap":"seismic","vskip":4})
# x,y temperature slice with velocity vectors
osiris.plot_slice(mydata.log_T,direction="z",vec=mydata.velocity,dx=100,axes=ax5,scalar_args={"cmap":"hot"},contour=mydata.level,contour_args={"fmt":"%i","label":True,"colors":"w","cmap":None,"levels":range(9,17)})

# Now update values with later snapshot
mydata.update_values(201)
# Re-plot x,y density slice with velocity vectors
osiris.plot_slice(mydata.log_rho,direction="z",vec=mydata.velocity,dx=100,axes=ax6)

fig.savefig("demo.pdf",bbox_inches="tight")
