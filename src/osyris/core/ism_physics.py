# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

"""
This file aims to re-introduce the ism_physics routines of osiris into Osyris.

To do:
-Opacities reader DONE
-Resistivities reader
-EOS reader DONE
"""

import struct
import os
import numpy as np
from ..core import Array
from .. import config
from .. import units
from ..io import utils
from scipy.interpolate import RegularGridInterpolator

def ism_interpolate(table_container=None, values=[0], points=[0], in_log=False):

	func = RegularGridInterpolator(table_container["grid"], values)

	if in_log:
		return func(points)
	else:
		return np.power(10.0, func(points))

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
	offsets = {"i":0, "n":0, "d":0}

	# Get table dimensions
	theTable["nx"] = np.array(utils.read_binary_data(fmt="3i",content=data,increment=False))

	# Read table coordinates:

	# x: density
	offsets["i"] += 3
	offsets["n"] += 9
	offsets["d"] += 1
	theTable["dens"] = utils.read_binary_data(fmt="%id"%theTable["nx"][0],content=data,offsets=offsets,increment=False)
	offsets["n"] -= 1

	# y: gas temperature
	offsets["n"] += theTable["nx"][0]
	offsets["d"] += 1
	theTable["tgas"] = utils.read_binary_data(fmt="%id"%theTable["nx"][1],content=data,offsets=offsets,increment=False)
	offsets["n"] -= 1

	# z: radiation temperature
	offsets["n"] += theTable["nx"][1]
	offsets["d"] += 1
	theTable["trad"] = utils.read_binary_data(fmt="%id"%theTable["nx"][2],content=data,offsets=offsets,increment=False)
	offsets["n"] -= 1

	# Now read opacities
	array_size = np.prod(theTable["nx"])
	array_fmt  = "%id" % array_size

	# Planck mean
	offsets["n"] += theTable["nx"][2]
	offsets["d"] += 1
	theTable["kappa_p"] = np.reshape(utils.read_binary_data(fmt=array_fmt,content=data, \
				offsets=offsets,increment=False),theTable["nx"],order="F")
	offsets["n"] -= 1

	# Rosseland mean
	offsets["n"] += array_size
	offsets["d"] += 1
	theTable["kappa_r"] = np.reshape(utils.read_binary_data(fmt=array_fmt,content=data, \
				offsets=offsets,increment=False),theTable["nx"],order="F")
	offsets["n"] -= 1

	del data

	theTable["grid"] = (theTable["dens"],theTable["tgas"],theTable["trad"])


	print("Opacity table read successfully")

	return theTable

def get_opacities(dataset, fname, variables={"kappa_p":"cm^2/g","kappa_r":"cm^2/g"}):
	"""
	Create opacity variables from interpolation of opacity table values in fname.
	"""

	if "opacity_table" not in dataset.meta:
		dataset.meta["opacity_table"] = read_opacity_table(fname=fname)

	if "radiative_temperature" not in dataset["hydro"]:
		print("Radiative temperature is not defined. Computing it now...", end="")
		dataset["hydro"]["radiative_temperature"] = values = (dataset["hydro"]["radiative_energy_1"]/units["radiation_constant"])**.25
		print(" done!")
	pts = np.array([np.log10(dataset["hydro"]["density"].values),np.log10(dataset["hydro"]["temperature"].values),np.log10(dataset["hydro"]["radiative_temperature"].values)]).T
	for var in variables:
		print("Interpolating "+var+"...", end="")
		vals = ism_interpolate(dataset.meta["opacity_table"], dataset.meta["opacity_table"][var], pts)
		print(" done!")
		dataset["hydro"][var] = Array(values = vals, unit = variables[var])

	return

def read_eos_table(fname):
	"""
	Read binary EOS table in fname
	"""

	print("Loading EOS table: "+fname)

	# Read binary EOS file
	with open(fname, mode='rb') as f:
	    data = f.read()

	# Define data fields. Note that the order is important!
	data_fields = ["rho_eos","ener_eos","temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]

	# Create table container
	theTable = dict()

	# Initialise offset counters and start reading data
	offsets = {"i":0, "n":0, "d":0}

	# Get table dimensions
	theTable["nx"] = np.array(utils.read_binary_data(fmt="2i",content=data, increment=False))

	# Get table limits
	offsets["i"] += 2
	offsets["d"] += 1
	[theTable["rhomin"],theTable["rhomax"],theTable["emin"],theTable["emax"],theTable["yHe"]] = \
		utils.read_binary_data(fmt="5d",content=data,offsets=offsets, increment=False)
	offsets["n"] -= 1

	array_size = np.prod(theTable["nx"])
	array_fmt  = "%id" % array_size
	offsets["n"] += 5
	offsets["d"] += 1

	# Now loop through all the data fields
	for i in range(len(data_fields)):
		theTable[data_fields[i]] = np.reshape(utils.read_binary_data(fmt=array_fmt,content=data, \
			offsets=offsets, increment=False),theTable["nx"],order="F")
		offsets["n"] += array_size
		offsets["n"] -= 1
		offsets["d"] += 1

	del data

	Eint = theTable["ener_eos"]/theTable["rho_eos"]
	theTable["grid"] = (np.log10(theTable["rho_eos"][:,0]), np.log10(Eint[0,:]))

	print("EOS table read successfully")

	return theTable

def get_eos(dataset, fname, variables={"temp_eos":"K","pres_eos":"dyn/cm^2","s_eos":"erg/K/g","cs_eos":"cm/s","xH_eos":None,"xH2_eos":None,"xHe_eos":None,"xHep_eos":None}):
	"""
	Create EOS variables from interpolation of eos table values in fname.
	"""

	if dataset.meta["eos"] == 0:
		print("Simulation data did not use a tabulated EOS. Exiting.")
		return
	if "eos_table" not in dataset.meta:
		dataset.meta["eos_table"] = read_eos_table(fname=fname)

	pts = np.array([np.log10(dataset["hydro"]["density"].values), np.log10(dataset["hydro"]["internal_energy"].values/dataset["hydro"]["density"].values)]).T
	for var in variables:
		print("Interpolating "+var+"...", end="")
		vals = ism_interpolate(dataset.meta["eos_table"],np.log10(dataset.meta["eos_table"][var]),pts)
		dataset["hydro"][var] = Array(values = vals, unit = variables[var])
		print(" done!")


