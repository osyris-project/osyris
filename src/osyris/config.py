# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np

#===================================================================================
# Define default values so that you don't have to specify them every time.
#===================================================================================
default_values = {
    "nout"        : 1    ,
    "lmax"        : 0    ,
    "center"      : None ,
    "dx"          : 0.0  ,
    "dy"          : 0.0  ,
    "dz"          : 0.0  ,
    "scale"       : "au" ,
    "time_unit"   : "kyr",
    "verbose"     : False,
    "path"        : ""   ,
    "variables"   : []   ,
    "colormap"    : "viridis"
}

#===================================================================================
# Define default units for some selected variables. The syntax is:
# "my_variable_name" : ["combination of unit_d, unit_l, unit_t" , "string to be displayed on grid axes"]
#===================================================================================
default_units = {
    ## Example for internal energy
    #"passive_scalar_4" : ["unit_d*((unit_l/unit_t)**2)" , "erg/cm3"],
}

#=======================================================================================
# Common variables
#=======================================================================================
constants = {
    "cm"  : 1.0,
    "au"  : 1.495980e+13,
    "pc"  : 3.085678e+18,
    "s"   : 1.0,
    "yr"  : 365.25*86400.0,
    "kyr" : 365.25*86400.0*1000.0,
    "msun": 1.9889e33,
    "a_r" : 7.56591469318689378e-015,
    "c"   : 2.9979250e+10
}

#=======================================================================================
# Custom colormaps
#=======================================================================================
cmaps = {
    "osyris"  : ["#2b3c4e","#249593","#db6a6c","#ffffff"],
    "osyris2" : ["#2b3c4e","#249593","#ffffff","#db6a6c","#9e4d4e"],
    "osyris3" : ["#3d3d6b","#2a7b9b","#00baad","#57c785","#add45c","#ffc300","#ff8d1a","#ff5733","#c70039","#900c3f","#511849"],
    "osyris4" : ["#000000","#ff5b00","#ffff00","#00ff00","#2bc184","#3d3d6b","#ffffff","#0000ff"],
    "osyris5" : ["#ff0000","#ffff00","#00ff00","#00ffff","#0000ff","#000000","#ffffff"]
}

for key in cmaps.keys():
    cmap = LinearSegmentedColormap.from_list(key, cmaps[key])
    cmap_r = LinearSegmentedColormap.from_list(key+"_r", cmaps[key][::-1])
    plt.register_cmap(cmap=cmap)
    plt.register_cmap(cmap=cmap_r)

#===================================================================================
# Here are some additional variables that are to be computed every time data is
# loaded.
#===================================================================================
def additional_variables(holder):

    # Velocity field (in case conservative variables are dumped)
    holder.new_field(name="velocity",operation="momentum/density",unit="cm/s",label="velocity",verbose=False)

    # Magnetic field
    holder.new_field(name="B",operation="0.5*(B_left+B_right)",unit="G",label="B",verbose=False)

    # Mass and radius
    holder.new_field(name="r",operation="np.sqrt(x**2 + y**2 + z**2)",unit=holder.info["scale"],label="Radius",verbose=False)
    holder.new_field(name="mass",operation="density*((dx*"+str(constants[holder.info["scale"]])+")**3)/"+str(constants["msun"]),unit="Msun",label="Mass",verbose=False)

    # Commonly used log quantities
    holder.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)",verbose=False)
    holder.new_field(name="log_T",operation="np.log10(temperature)",unit="K",label="log(T)",verbose=False)
    holder.new_field(name="log_B",operation="np.log10(B.values)",unit="G",label="log(B)",verbose=False)
    holder.new_field(name="log_m",operation="np.log10(mass)",unit="Msun",label="log(Mass)",verbose=False)
    with np.errstate(divide="ignore"):
        holder.new_field(name="log_r",operation="np.log10(r)",unit=holder.info["scale"],label="log(Radius)",verbose=False)

    # Photon density (in case conservative variables are not dumped)
    if hasattr(holder,"photon_flux_density_1"):
        for igrp in range(holder.info_rt["nGroups"]):
            holder.new_field(name="photon_density_"+str(igrp+1),operation="photon_flux_density_"+str(igrp+1)+"/("+str(constants["c"]*holder.info_rt["rt_c_frac"])+")",unit="photons/cm3",label="photon density "+str(igrp+1),verbose=False)

    #========================== ADD YOUR VARIABLES HERE ============================

    return
