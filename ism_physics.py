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

from load_ramses_data import get_binary_data
import numpy as np
import struct
from scipy.interpolate import RegularGridInterpolator

#===================================================================================
# An empty EOS table class which will be filled up as the EOS file is read
#===================================================================================
class EosTable():
    
    def __init__(self):
                
        return

#===================================================================================
# Function to read in binary file containing EOS table
#===================================================================================
def read_eos_table(fname="tab_eos.dat"):
    
    print("Loading EOS table: "+fname)
    
    # Read binary EOS file
    with open(fname, mode='rb') as eos_file:
        eosContent = eos_file.read()
    eos_file.close()
    
    # Define data fields. Note that the order is important!
    data_fields = ["rho_eos","ener_eos","temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]

    # Create table container
    theTable = EosTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    # Get table dimensions
    [theTable.nRho,theTable.nEner] = get_binary_data(fmt="2i",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat,nstrin=nstrin,nquadr=nquadr)
    
    # Get table limits
    ninteg += 2
    nlines += 1
    [theTable.rhomin,theTable.rhomax,theTable.emin,theTable.emax,theTable.yHe] = \
        get_binary_data(fmt="5d",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat,nstrin=nstrin,nquadr=nquadr)
        
    array_size = theTable.nRho*theTable.nEner
    array_fmt  = "%id" % array_size
    nfloat += 5
    nlines += 1
    
    # Now loop through all the data fields
    for i in range(len(data_fields)):
        setattr(theTable,data_fields[i],np.reshape(get_binary_data(fmt=array_fmt,content=eosContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat,nstrin=nstrin,nquadr=nquadr),[theTable.nRho,theTable.nEner],order="F"))
        nfloat += array_size
        nlines += 1
    
    del eosContent
    
    print("EOS table read successfully")

    return theTable


#===================================================================================
# Function to read in binary file containing EOS table
#===================================================================================
def get_eos_variables(holder,eos_fname="tab_eos.dat",variables=["temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]):
    
    if holder.info["eos"] == 0:
        print("Simulation data did not use a tabulated EOS. Exiting")
        return
    else:
        try:
            n = holder.eos_table.nRho
        except AttributeError:
            holder.eos_table = read_eos_table(fname=eos_fname)

        for i in range(len(variables)):
            print("Interpolating "+variables[i])
            grid = (np.log10(holder.eos_table.rho_eos[:,0]), np.log10(holder.eos_table.ener_eos[0,:]/holder.eos_table.rho_eos[0,:]))
            func = RegularGridInterpolator(grid, np.log10(getattr(holder.eos_table,variables[i])))
            pts  = np.array([np.log10(holder.density.values), np.log10(holder.internal_energy.values)]).T
            holder.new_field(name=variables[i],label=variables[i],values=np.power(10.0,func(pts)),verbose=False,norm=1.0)

    return
