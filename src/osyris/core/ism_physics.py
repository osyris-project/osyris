# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

"""
This file aims to re-introduce the ism_physics routines of osiris into Osyris.

To do:
-Opacities reader
-Resistivities reader
-EOS reader
"""

import os
import numpy as np
from ..core import Array
from .. import config
from .. import units
from ..io import utils
from scipy.interpolate import RegularGridInterpolator

def read_opacity_table(fname):
	"""
	Read binary opacity table in fname.
	"""

	print("Loading opacity table: "+fname)

	with open(fname, "rb") as f:
		data = f.read()

	# Create table container
	theTable = dict()

	# Initialise offset counters and start reading data
	offsets = {"i":0, "f":0, "l":0}

	# Get table dimensions
	theTable["nx"] = np.array(utils.read_binary_data(fmt="3i",content=data))

	# Read table coordinates:

	# x: density
	offsets["i"] += 3
	offsets["n"] += 9
	offsets["d"] += 1
	theTable["dens"] = utils.read_binary_data(fmt="%id"%theTable["nx"][0],content=data,offsets=offsets)

	# y: gas temperature
	offsets["n"] += theTable["nx"][0]
	offsets["d"] += 1
	theTable["tgas"] = utils.read_binary_data(fmt="%id"%theTable["nx"][1],content=data,offsets=offsets)

	# z: radiation temperature
	offsets["n"] += theTable["nx"][1]
	offsets["d"] += 1
	theTable["trad"] = utils.read_binary_data(fmt="%id"%theTable["nx"][2],content=data,offsets=offsets)

	# Now read opacities
	array_size = np.prod(theTable["nx"])
	array_fmt  = "%id" % array_size

	#print theTable.nx,theTable.ny,theTable.nz

	# Planck mean
	offsets["n"] += theTable["nx"][2]
	offsets["d"] += 1
	theTable.kappa_p = np.reshape(utils.read_binary_data(fmt=array_fmt,content=data, \
	            offsets=offsets),theTable["nx"],order="F")

	# Rosseland mean
	offsets["n"] += array_size
	offsets["d"] += 1
	theTable["kappa_r"] = np.reshape(utils.read_binary_data(fmt=array_fmt,content=data, \
	            offsets=offsets),theTable["nx"],order="F")

	del data

	theTable["grid"] = (theTable["dens"],theTable["tgas"],theTable["trad"])


	print("Opacity table read successfully")

	return theTable
"""
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
"""