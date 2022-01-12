# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

"""
This file aims to re-introduce the ism_physics routines of osiris into Osyris.

To do:
-Opacities reader DONE
-Resistivities reader
-EOS reader
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

def read_binary_data(fmt="", offsets=None, content=None, correction=0):

	if offsets is not None:
		ninteg = offsets["i"]
		nfloat = offsets["n"]
		nlines = offsets["d"]
	else:
		ninteg = 0
		nfloat = 0
		nlines = 0
	nstrin = 0
	nquadr = 0
	nlongi = 0

	offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16 + 4 + correction
	byte_size = {"i":4,"d":8,"q":8}
	if len(fmt) == 1:
	    mult = 1
	else:
	    mult = eval(fmt[0:len(fmt)-1])
	pack_size = mult*byte_size[fmt[-1]]

	return struct.unpack(fmt, content[offset:offset+pack_size])


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
	theTable["nx"] = np.array(read_binary_data(fmt="3i",content=data))

	# Read table coordinates:

	# x: density
	offsets["i"] += 3
	offsets["n"] += 9
	offsets["d"] += 1
	theTable["dens"] = read_binary_data(fmt="%id"%theTable["nx"][0],content=data,offsets=offsets)

	# y: gas temperature
	offsets["n"] += theTable["nx"][0]
	offsets["d"] += 1
	theTable["tgas"] = read_binary_data(fmt="%id"%theTable["nx"][1],content=data,offsets=offsets)

	# z: radiation temperature
	offsets["n"] += theTable["nx"][1]
	offsets["d"] += 1
	theTable["trad"] = read_binary_data(fmt="%id"%theTable["nx"][2],content=data,offsets=offsets)

	# Now read opacities
	array_size = np.prod(theTable["nx"])
	array_fmt  = "%id" % array_size

	#print theTable.nx,theTable.ny,theTable.nz

	# Planck mean
	offsets["n"] += theTable["nx"][2]
	offsets["d"] += 1
	theTable["kappa_p"] = np.reshape(read_binary_data(fmt=array_fmt,content=data, \
	            offsets=offsets),theTable["nx"],order="F")

	# Rosseland mean
	offsets["n"] += array_size
	offsets["d"] += 1
	theTable["kappa_r"] = np.reshape(read_binary_data(fmt=array_fmt,content=data, \
	            offsets=offsets),theTable["nx"],order="F")

	del data

	theTable["grid"] = (theTable["dens"],theTable["tgas"],theTable["trad"])


	print("Opacity table read successfully")

	return theTable

def get_opacities(dataset, fname, variables=["kappa_p","kappa_r"]):

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
		dataset["hydro"][var] = Array(values = vals, unit = "cm*cm/g")

	return
