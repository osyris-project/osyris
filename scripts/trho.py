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

import osiris
import sys

fname = ""
if len(sys.argv) > 2:
    fname = sys.argv[2]
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = -1

# Load data
mydata = osiris.RamsesData(nout=nout,center="max:density",scale="au")

# Density vs T
if(len(fname) > 0):
    osiris.plot_histogram(mydata.log_rho,mydata.log_T,fname=fname,scalar_args={"cmap":"osiris_r,log"},outline=True)
else:
    osiris.plot_histogram(mydata.log_rho,mydata.log_T,block=True,scalar_args={"cmap":"osiris_r,log"},outline=True)
