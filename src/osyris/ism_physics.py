# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

from .engine import get_binary_data
from . import config as conf
import numpy as np
from scipy.interpolate import RegularGridInterpolator

#===================================================================================
# An empty container class which will be filled up as the EOS, opacity, resistivity
# tables are read.
#===================================================================================
class IsmTable():

    def __init__(self):

        return


#===================================================================================
# Generic function to interpolate simulation data onto opacity table
#===================================================================================
def ism_interpolate(table_container=None,values=[0],points=[0],in_log=False):

    func = RegularGridInterpolator(table_container.grid,values)

    if in_log:
        return func(points)
    else:
        return np.power(10.0,func(points))


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
    theTable.nx = np.array(get_binary_data(fmt="2i",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat))

    # Get table limits
    ninteg += 2
    nlines += 1
    [theTable.rhomin,theTable.rhomax,theTable.emin,theTable.emax,theTable.yHe] = \
        get_binary_data(fmt="5d",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    array_size = np.prod(theTable.nx)
    array_fmt  = "%id" % array_size
    nfloat += 5
    nlines += 1

    # Now loop through all the data fields
    for i in range(len(data_fields)):
        setattr(theTable,data_fields[i],np.reshape(get_binary_data(fmt=array_fmt,content=eosContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F"))
        nfloat += array_size
        nlines += 1

    del eosContent

    #theTable.grid = (np.log10(theTable.rho_eos[:,0]), np.log10(theTable.ener_eos[0,:]/theTable.rho_eos[0,:]))

    Eint = theTable.ener_eos/theTable.rho_eos
    theTable.grid = (np.log10(theTable.rho_eos[:,0]), np.log10(Eint[0,:]))

    #print theTable.rho_eos[:,0]
    #print Eint[0,:]
    #print theTable.ener_eos[0,:]

    print("EOS table read successfully")

    return theTable


#===================================================================================
# Function to interpolate simulation data onto EOS table
#===================================================================================
def get_eos(holder,fname="tab_eos.dat",variables=["temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]):

    if holder.info["eos"] == 0:
        print("Simulation data did not use a tabulated EOS. Exiting")
        return
    else:
        try:
            n = holder.eos_table.nRho
        except AttributeError:
            holder.eos_table = read_eos_table(fname=fname)

        pts = np.array([np.log10(holder.density.values), np.log10(holder.internal_energy.values/holder.density.values)]).T
        for var in variables:
            print("Interpolating "+var)
            vals = ism_interpolate(holder.eos_table,np.log10(getattr(holder.eos_table,var)),pts)
            holder.new_field(name=var,label=var,values=vals,verbose=False)

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
    theTable.nx = np.array(get_binary_data(fmt="3i",content=kappaContent))

    ## Get table limits
    #[theTable.dx,theTable.dy,theTable.dz,theTable.xmin,theTable.xmax,theTable.ymin,theTable.ymax,theTable.zmin,theTable.zmax] = \
        #get_binary_data(fmt="9d",content=kappaContent,correction=12)

    # Read table coordinates:

    # x: density
    ninteg += 3
    nfloat += 9
    nlines += 1
    theTable.dens = get_binary_data(fmt="%id"%theTable.nx[0],content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    # y: gas temperature
    nfloat += theTable.nx[0]
    nlines += 1
    theTable.tgas = get_binary_data(fmt="%id"%theTable.nx[1],content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    # z: radiation temperature
    nfloat += theTable.nx[1]
    nlines += 1
    theTable.trad = get_binary_data(fmt="%id"%theTable.nx[2],content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    # Now read opacities
    array_size = np.prod(theTable.nx)
    array_fmt  = "%id" % array_size

    #print theTable.nx,theTable.ny,theTable.nz

    # Planck mean
    nfloat += theTable.nx[2]
    nlines += 1
    theTable.kappa_p = np.reshape(get_binary_data(fmt=array_fmt,content=kappaContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")

    # Rosseland mean
    nfloat += array_size
    nlines += 1
    theTable.kappa_r = np.reshape(get_binary_data(fmt=array_fmt,content=kappaContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")

    del kappaContent

    theTable.grid = (theTable.dens,theTable.tgas,theTable.trad)

    #print theTable.dens
    #print np.shape(theTable.dens)
    ##print theTable.tgas
    ##print theTable.trad


    print("Opacity table read successfully")

    return theTable


#===================================================================================
# Function to interpolate simulation data onto opacity table
#===================================================================================
def get_opacities(holder,fname="vaytet_grey_opacities3D.bin",variables=["kappa_p","kappa_r"]):

    try:
        n = holder.opacity_table.nx
    except AttributeError:
        holder.opacity_table = read_opacity_table(fname=fname)

    if not hasattr(holder,"radiative_temperature"):
        print("Radiative temperature is not defined. Computing it now.")
        holder.new_field(name="radiative_temperature",operation="(radiative_energy_1/"+str(conf.constants["a_r"])+")**0.25",label="Trad")

    pts = np.array([np.log10(holder.density.values),np.log10(holder.temperature.values),np.log10(holder.radiative_temperature.values)]).T
    for var in variables:
        print("Interpolating "+var)
        vals = ism_interpolate(holder.opacity_table,getattr(holder.opacity_table,var),pts)
        holder.new_field(name=var,label=var,values=vals,verbose=False,unit="cm2/g")

    return


#===================================================================================
#===================================================================================
#===================================================================================

#===================================================================================
# Function to read in binary file containing resistivity table
#===================================================================================
def read_resistivity_table(fname="resistivities_masson2016.bin"):

    print("Loading resistivity table: "+fname)

    # Read binary resistivity file
    with open(fname, mode='rb') as res_file:
        resContent = res_file.read()
    res_file.close()

    # Create table container
    theTable = IsmTable()

    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0

    # Get length of record on first line to determine number of dimensions in table
    rec_size = get_binary_data(fmt="i",content=resContent,correction=-4)
    theTable.ndims = int(rec_size[0]/4)

    # Get table dimensions
    theTable.nx = np.array(get_binary_data(fmt="%ii"%theTable.ndims,content=resContent))

    # Read table coordinates:

    # 1: density
    ninteg += theTable.ndims
    nlines += 1
    theTable.dens = get_binary_data(fmt="%id"%theTable.nx[0],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    # 2: gas temperature
    nfloat += theTable.nx[0]
    nlines += 1
    theTable.tgas = get_binary_data(fmt="%id"%theTable.nx[1],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    if theTable.ndims == 4:
        # 3: ionisation rate
        nfloat += theTable.nx[1]
        nlines += 1
        theTable.ionx = get_binary_data(fmt="%id"%theTable.nx[2],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    # 4: magnetic field
    nfloat += theTable.nx[-2]
    nlines += 1
    theTable.bmag = get_binary_data(fmt="%id"%theTable.nx[-1],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

    # Now read resistivities
    array_size = np.prod(theTable.nx)
    array_fmt  = "%id" % array_size

    # Ohmic resistivity
    nfloat += theTable.nx[-1]
    nlines += 1
    theTable.eta_ohm = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")

    # Ambipolar resistivity
    nfloat += array_size
    nlines += 1
    theTable.eta_ad = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")

    # Hall resistivity
    nfloat += array_size
    nlines += 1
    theTable.eta_hall = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")

    # Hall sign
    nfloat += array_size
    nlines += 1
    theTable.eta_hsig = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")

    del resContent

    if theTable.ndims == 4:
        theTable.grid = (theTable.dens,theTable.tgas,theTable.ionx,theTable.bmag)
    elif theTable.ndims == 3:
        theTable.grid = (theTable.dens,theTable.tgas,theTable.bmag)

    # Additional parameters
    theTable.scale_dens = 0.844*2.0/1.66e-24 # 2.0*H2_fraction/mH
    theTable.ionis_rate = 1.0e-17

    print("Resistivity table read successfully")

    return theTable


#===================================================================================
# Function to interpolate simulation data onto resistivity table
#===================================================================================
def get_resistivities(holder,fname="resistivities_masson2016.bin",variables=["eta_ohm","eta_ad","eta_hall"]):

    try:
        n = holder.resistivity_table.nx
    except AttributeError:
        holder.resistivity_table = read_resistivity_table(fname=fname)

    try:
        rho_to_nH = holder.resistivity_table.scale_dens/holder.info["mu_gas"]
    except KeyError:
        rho_to_nH = holder.resistivity_table.scale_dens/2.31 # Default mean atomic weight

    if holder.resistivity_table.ndims == 4:
        if not hasattr(holder,"ionisation_rate"):
            holder.new_field(name="ionisation_rate",values=[holder.resistivity_table.ionis_rate]*holder.info["ncells"],label="Ionisation rate")
        pts = np.array([np.log10(holder.density.values*rho_to_nH),np.log10(holder.temperature.values), \
                        np.log10(holder.ionisation_rate.values),np.log10(holder.B.values)]).T
    elif holder.resistivity_table.ndims == 3:
        pts = np.array([np.log10(holder.density.values*rho_to_nH),np.log10(holder.temperature.values), \
                        np.log10(holder.B.values)]).T

    for var in variables:
        print("Interpolating "+var)
        vals = ism_interpolate(holder.resistivity_table,getattr(holder.resistivity_table,var),pts)
        if var == "eta_hall":
            hall_sign = np.sign(ism_interpolate(holder.resistivity_table,holder.resistivity_table.eta_hsig,pts,in_log=True))
            holder.new_field(name=var,label=var,values=vals*hall_sign,verbose=False,unit="s")
        else:
            holder.new_field(name=var,label=var,values=vals,verbose=False,unit="s")

    return
