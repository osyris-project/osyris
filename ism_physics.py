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
# An empty container class which will be filled up as the EOS, opacity, resistivity
# tables are read.
#===================================================================================
class IsmTable():
    
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
    theTable = IsmTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    # Get table dimensions
    [theTable.nRho,theTable.nEner] = get_binary_data(fmt="2i",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # Get table limits
    ninteg += 2
    nlines += 1
    [theTable.rhomin,theTable.rhomax,theTable.emin,theTable.emax,theTable.yHe] = \
        get_binary_data(fmt="5d",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
        
    array_size = theTable.nRho*theTable.nEner
    array_fmt  = "%id" % array_size
    nfloat += 5
    nlines += 1
    
    # Now loop through all the data fields
    for i in range(len(data_fields)):
        setattr(theTable,data_fields[i],np.reshape(get_binary_data(fmt=array_fmt,content=eosContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),[theTable.nRho,theTable.nEner],order="F"))
        nfloat += array_size
        nlines += 1
    
    del eosContent
    
    print("EOS table read successfully")

    return theTable


#===================================================================================
# Function to interpolate simulation data onto EOS table
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







#===================================================================================
#===================================================================================
#===================================================================================

#===================================================================================
# Function to read in binary file containing opacity table
#===================================================================================
def read_opacity_table(fname="vaytet_grey_opacities3D.bin"):
    
    print("Loading opacity table: "+fname)
    
    # Read binary resistivity file
    with open(fname, mode='rb') as kappa_file:
        kappaContent = kappa_file.read()
    kappa_file.close()
    
    # Create table container
    theTable = IsmTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    # Get table dimensions
    [theTable.nx,theTable.ny,theTable.nz] = get_binary_data(fmt="3i",content=kappaContent)
    
    # Get table limits
    [theTable.dx,theTable.dy,theTable.dz,theTable.xmin,theTable.xmax,theTable.ymin,theTable.ymax,theTable.zmin,theTable.zmax] = \
        get_binary_data(fmt="9d",content=kappaContent,correction=12)
    
    # Read table coordinates:
    
    # x: density
    ninteg += 3
    nfloat += 9
    nlines += 1
    theTable.dens = get_binary_data(fmt="%id"%theTable.nx,content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # y: gas temperature
    nfloat += theTable.nx
    nlines += 1
    theTable.tgas = get_binary_data(fmt="%id"%theTable.ny,content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # z: radiation temperature
    nfloat += theTable.ny
    nlines += 1
    theTable.trad = get_binary_data(fmt="%id"%theTable.nz,content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # Now read opacities
    array_size = theTable.nx*theTable.ny*theTable.nz
    array_fmt  = "%id" % array_size
    
    #print theTable.nx,theTable.ny,theTable.nz
    
    # Planck mean
    nfloat += theTable.nz
    nlines += 1
    theTable.kappa_p = np.reshape(get_binary_data(fmt=array_fmt,content=kappaContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),[theTable.nx,theTable.ny,theTable.nz],order="F")
    
    # Rosseland mean
    nfloat += array_size
    nlines += 1
    theTable.kappa_r = np.reshape(get_binary_data(fmt=array_fmt,content=kappaContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),[theTable.nx,theTable.ny,theTable.nz],order="F")
    
    del kappaContent
    print("Opacity table read successfully")
    
    return theTable


#===================================================================================
# Function to interpolate simulation data onto opacity table
#===================================================================================
def get_opacities(holder,opacity_fname="vaytet_grey_opacities3D.bin",variables=["kappa_p","kappa_r"]):
    
    
    try:
        n = holder.opacity_table.nx
    except AttributeError:
        holder.opacity_table = read_opacity_table(fname=opacity_fname)

    for i in range(len(variables)):
        print("Interpolating "+variables[i])
        grid = (holder.opacity_table.dens,holder.opacity_table.tgas,holder.opacity_table.trad)
        func = RegularGridInterpolator(grid, getattr(holder.opacity_table,variables[i]))
        #print 'got to here'
        #print np.log10(holder.radiative_energy_1.values)
        pts  = np.array([np.log10(holder.density.values),np.log10(holder.temperature.values),np.log10(holder.radiative_temperature.values)]).T
        holder.new_field(name=variables[i],label=variables[i],values=np.power(10.0,func(pts)),verbose=False,norm=1.0)

    return



#===================================================================================
#===================================================================================
#===================================================================================

#===================================================================================
# Function to read in binary file containing resistivity table
#===================================================================================
def read_resistivity_table(fname="resnh.bin"):
    
    print("Loading resistivity table: "+fname)
    
    # Read binary resistivity file
    with open(fname, mode='rb') as eta_file:
        etaContent = eta_file.read()
    eta_file.close()
    
    # Define data fields. Note that the order is important!
    #data_fields = ["rho_eos","ener_eos","temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]

    # Create table container
    theTable = IsmTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    #nx = struct.unpack("i", etaContent[0:4])
    #print nx
    
    # Get length of record on first line to determine number of dimensions in table
    rec_size = get_binary_data(fmt="i",content=etaContent,correction=-4)
    ndims = rec_size[0]/4
    
    # Get table dimensions
    theTable.nx = np.roll(np.array(get_binary_data(fmt="%ii"%ndims,content=etaContent)),1)
    print theTable.nx
    
    if ndims == 3:
        theTable.nx[0] += 2
    elif ndims == 4:
        theTable.nx[0] += 7
    
    # Now read the bulk of the table in one go
    ninteg += ndims
    nlines += 1
    theTable.res_chimie = np.reshape(get_binary_data(fmt="%id"%(np.prod(theTable.nx)),content=etaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")
    
    print theTable.res_chimie[:,0,0]
    print np.shape(theTable.res_chimie)
    

    return theTable















#open(42,file='resnh.dat', status='old')
        #read(42,*) nchimie, tchimie, nvarchimie
        #read(42,*)
        #read(42,*)
        #allocate(resistivite_chimie_x(-1:nvarchimie,nchimie,tchimie,1))
        #do i=1,tchimie
           #do j=1,nchimie
              #read(42,*)resistivite_chimie_x(0:nvarchimie,j,i,1),dummy,dummy,dummy,dummy,resistivite_chimie_x(-1,j,i,1)
#!              print *, resistivite_chimie_x(:,j,i)
           #end do
           #read(42,*)
        #end do
        #close(42)
        #rho_threshold=max(rho_threshold,resistivite_chimie_x(0,1,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
        #nminchimie=(resistivite_chimie_x(0,1,1,1))
        #dnchimie=(log10(resistivite_chimie_x(0,nchimie,1,1))-log10(resistivite_chimie_x(0,1,1,1)))/&
                 #&(nchimie-1)
#!                 print*, dnchimie,15.d0/50.d0
        #tminchimie=(resistivite_chimie_x(-1,1,1,1))
        #dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie,1))-log10(resistivite_chimie_x(-1,1,1,1)))/&
                 #&(tchimie-1)
#!                 print*, dtchimie,3.d0/50.d0
#!         close(333)
        #call rq
        #call nimhd_3dtable
     #else if(use_x3d==1)then

        #open(42,file='marchand2016_table.dat',form='unformatted')
        #read(42) nchimie, tchimie, xichimie, nvarchimie
        #allocate(resistivite_chimie_x(-2:nvarchimie+4,nchimie,tchimie,xichimie))
        #read(42) resistivite_chimie_x
        #close(42)

        #rho_threshold=max(rho_threshold,resistivite_chimie_x(-2,1,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
        #nminchimie=(resistivite_chimie_x(-2,1,1,1))
        #dnchimie=(log10(resistivite_chimie_x(-2,nchimie,1,1))-log10(resistivite_chimie_x(-2,1,1,1)))/&
                 #&(nchimie-1)
#!                 print*, dnchimie,15.d0/50.d0
        #tminchimie=(resistivite_chimie_x(-1,1,1,1))
        #dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie,1))-log10(resistivite_chimie_x(-1,1,1,1)))/&
                 #&(tchimie-1)
#!                 print*, dtchimie,3.d0/50.d0
        #ximinchimie=(resistivite_chimie_x(0,1,1,1))
        #dxichimie=(log10(resistivite_chimie_x(0,1,1,xichimie))-log10(resistivite_chimie_x(0,1,1,1)))/&
                 #&(xichimie-1)
        #call rq_3d
        #call nimhd_4dtable
     #else
        #print*, 'must choose an input for abundances or resistivities'
        #stop
     #endif
