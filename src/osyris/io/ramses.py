# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
import pandas as pd
import struct
import glob
from .. import config as conf
from .. import engine as eng

divider = "============================================"

def generate_fname(nout,path="",ftype="",cpuid=1,ext=""):

    if len(path) > 0:
        if path[-1] != "/":
            path=path+"/"

    if nout == -1:
        filelist = sorted(glob.glob(path+"output*"))
        number = filelist[-1].split("_")[-1]
    else:
        number = str(nout).zfill(5)

    infile = path+"output_"+number
    if len(ftype) > 0:
        infile += "/"+ftype+"_"+number
        if cpuid >= 0:
            infile += ".out"+str(cpuid).zfill(5)

    if len(ext) > 0:
        infile += ext

    return infile


# def read_parameter_file(fname="",dict_name="",evaluate=True,verbose=False,delimiter="="):

#     # Read info file and create dictionary
#     try:
#         with open(fname) as f:
#             content = f.readlines()
#         f.close()
#     except IOError:
#         # Clean exit if the file was not found
#         if verbose:
#             print("File not found: "+fname)
#         #if raise_error:
#             #raise IOError
#         #else:
#         return 0

#     setattr(self,dict_name,dict())
#     for line in content:
#         sp = line.split(delimiter)
#         if len(sp) > 1:
#             if evaluate:
#                 try:
#                     getattr(self,dict_name)[sp[0].strip()] = eval(sp[1].strip())
#                 except (NameError,SyntaxError):
#                     getattr(self,dict_name)[sp[0].strip()] = sp[1].strip()
#             else:
#                 getattr(self,dict_name)[sp[0].strip()] = sp[1].strip()

#     return 1



def load(nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="cm",path="",
                    update=False,variables=[]):


    df = pd.DataFrame()

    # Generate directory name from output number
    infile = generate_fname(nout,path)

    # Read info file and create info dictionary
    infofile = infile+"/info_"+infile.split("_")[-1]+".txt"

    info = pd.read_csv(infofile, delimiter="=", header=None)
    for key, value in zip(info[0], info[1]):
        try:
            x = float(value)
            if np.isfinite(x):
                df.attrs[key.strip()] = x
        except ValueError:
            df.attrs[key.strip()] = value

    print(df.attrs)

    # Add additional information
    df.attrs["center"   ] = center
    df.attrs["scale"    ] = scale
    df.attrs["infile"   ] = infile
    df.attrs["path"     ] = path
    df.attrs["boxsize"  ] = df.attrs["boxlen"]*df.attrs["unit_l"]
    df.attrs["time"     ] = df.attrs["time"]*df.attrs["unit_t"]
    df.attrs["dx_load"  ] = dx
    df.attrs["dy_load"  ] = dy
    df.attrs["dz_load"  ] = dz
    df.attrs["lmax"     ] = lmax
    df.attrs["variables"] = variables

    # Convert to integers
    df.attrs["ncpu"]         = int(df.attrs["ncpu"]        )
    df.attrs["ndim"]         = int(df.attrs["ndim"]        )
    df.attrs["levelmin"]     = int(df.attrs["levelmin"]    )
    df.attrs["levelmax"]     = int(df.attrs["levelmax"]    )
    df.attrs["ngridmax"]     = int(df.attrs["ngridmax"]    )
    df.attrs["nstep_coarse"] = int(df.attrs["nstep_coarse"])
    df.attrs["ngrp"]         = int(df.attrs["ngrp"]        )
    df.attrs["ir_cloud"]     = int(df.attrs["ir_cloud"]        )
    df.attrs["eos"]          = int(df.attrs["eos"]        )

       # if nout==-1:
    #     df.attrs["nout" ] = int(infile.split("_")[-1])
    # else:
    #     df.attrs["nout" ] = nout

    # # Read namelist file and create namelist dictionary
    # nmlfile = infile+"/namelist.txt"
    # self.read_parameter_file(fname=nmlfile,dict_name="namelist",evaluate=False)

    # print(divider)

    # Now go through all the variables and check if they are to be read or skipped
    list_vars   = []
    var_read    = []
    var_group   = []
    var_type    = []
    xyz_strings = "xyz"

    # Start with hydro variables ==================================

    # Read the number of variables from the hydro_file_descriptor.txt
    # and select the ones to be read if specified by user
    hydrofile = infile+"/hydro_file_descriptor.txt"

    descriptor = pd.read_csv(hydrofile, header=None, comment="#")

    var_read = [True] * len(descriptor)
    list_vars = list(descriptor[1])
    var_type = list(descriptor[2])

    # var_read.append(True)
    #             list_vars.append(v)
    #             var_group.append("hydro")
    #             var_type.append(t)



    hydro = True
    df.attrs["nvar_hydro"] = len(descriptor)
    # try:
    #     with open(hydrofile) as f:
    #         content = f.readlines()
    # except IOError:
    #     hydro = False
    # if hydro:
    #     # Store the total number of hydro variables
    #     df.attrs["nvar_hydro"] = len(content) - 2
    #     # Now add to the list of variables to be read
    #     for line in content[2:]:
    #         sp = line.split(",")
    #         v = sp[1].strip()
    #         t = sp[2].strip()
    #         if (len(variables) == 0) or (v in variables) or ("hydro" in variables):
    #             var_read.append(True)
    #             list_vars.append(v)
    #             var_group.append("hydro")
    #             var_type.append(t)
    #         else:
    #             var_read.append(False)

    # # Now for gravity ==================================

    # # Check if self-gravity files exist
    # grav_fname = self.generate_fname(nout,path,ftype="grav",cpuid=1)
    gravity = False
    # df.attrs["nvar_grav"] = 0
    # try:
    #     with open(grav_fname, mode='rb') as grav_file:
    #         gravContent = grav_file.read()
    # except IOError:
    #     gravity = False

    # # Add gravity fields
    # if gravity:
    #     df.attrs["nvar_grav"] = 4
    #     content = ["grav_potential"]
    #     for n in range(df.attrs["ndim"]):
    #         content.append("grav_acceleration_"+xyz_strings[n])

    #     # Now add to the list of variables to be read
    #     for line in content:
    #         if (len(variables) == 0) or (line.strip() in variables) or ("gravity" in variables) or ("grav" in variables):
    #             var_read.append(True)
    #             list_vars.append(line.strip())
    #             var_group.append("grav")
    #         else:
    #             var_read.append(False)

    # # Now for rt ==================================

    # rtfile = infile+"/rt_file_descriptor.txt"
    rt = False
    # df.attrs["nvar_rt"] = 0
    # try:
    #     with open(rtfile) as f:
    #         content = f.readlines()
    #     # f.close()
    # except IOError:
    #     rt = False
    # if rt:
    #     # Store the total number of rt variables
    #     df.attrs["nvar_rt"] = len(content) - 2
    #     # Now add to the list of variables to be read
    #     for line in content[2:]:
    #         sp = line.split(",")
    #         v = sp[1].strip()
    #         t = sp[2].strip()
    #         if (len(variables) == 0) or (v in variables) or ("rt" in variables):
    #             var_read.append(True)
    #             list_vars.append(v)
    #             var_group.append("rt")
    #             var_type.append(t)
    #         else:
    #             var_read.append(False)

    # Make sure we always read the coordinates
    var_amr = ["level","x","y","z","dx","cpu"]
    list_vars.extend(var_amr)
    var_read.extend([True] * len(var_amr))
    # var_group.extend(["amr"] * len(var_amr))
    nvar_read = len(list_vars)
    print(nvar_read)

    # # Now for particles ==================================

    # # TODO: refactor this code to use the same function for reading hydro,
    # # rt, and part file descriptors

    # particles = True
    # partfile = infile+"/part_file_descriptor.txt"
    # df.attrs["npart_tot"] = 0
    # try:
    #     with open(partfile) as f:
    #         content = f.readlines()
    #     # f.close()
    # except IOError:
    #     particles = False
    # if particles:
    #     part_read = []
    #     part_vars = []
    #     part_type = []
    #     # Store the total number of part variables
    #     df.attrs["nvar_part"] = len(content) - 2
    #     # Now add to the list of variables to be read
    #     for line in content[2:]:
    #         sp = line.split(",")
    #         v = "part_" + sp[1].strip()
    #         t = sp[2].strip()
    #         if (len(variables) == 0) or (v in variables) or ("part" in variables):
    #             part_read.append(True)
    #             part_vars.append(v)
    #             part_type.append(t)
    #         else:
    #             part_read.append(False)

    # # Load sink particles if any
    # self.read_sinks()

    # Find the center
    # xc,yc,zc = self.find_center(dx,dy,dz)
    xc,yc,zc = 0.5, 0.5, 0.5

    # Now read the amr and hydro files =============================================
    # We have to open the files in binary format, and count all the bytes in the ===
    # file structure to extract just the data we need. =============================
    # See output_amr.f90 and output_hydro.f90 in the RAMSES source. ================
    print("Processing %i files in " % (df.attrs["ncpu"]) + infile)

    # Define the size of the region to be read
    lconvert = conf.constants[scale]/(df.attrs["boxlen"]*df.attrs["unit_l"])
    if dx > 0.0:
        xmin = xc - 0.5*dx*lconvert
        xmax = xc + 0.5*dx*lconvert
    else:
        xmin = 0.0
        xmax = 1.0
    if dy > 0.0:
        ymin = yc - 0.5*dy*lconvert
        ymax = yc + 0.5*dy*lconvert
    else:
        ymin = 0.0
        ymax = 1.0
    if dz > 0.0:
        zmin = zc - 0.5*dz*lconvert
        zmax = zc + 0.5*dz*lconvert
    else:
        zmin = 0.0
        zmax = 1.0

    if lmax==0:
       lmax = df.attrs["levelmax"]

    # We will store the cells in a dictionary which we build as we go along.
    # The final concatenation into a single array will be done once at the end.
    data_pieces = dict()
    npieces = 0
    part_pieces = dict()
    npieces_part = 0
    npart_count = 0

    # Allocate work arrays
    twotondim = 2**df.attrs["ndim"]
    xcent = np.zeros([8,3],dtype=np.float64)
    xg    = np.zeros([df.attrs["ngridmax"],3],dtype=np.float64)
    son   = np.zeros([df.attrs["ngridmax"],twotondim],dtype=np.int32)
    var   = np.zeros([df.attrs["ngridmax"],twotondim,nvar_read],dtype=np.float64)
    xyz   = np.zeros([df.attrs["ngridmax"],twotondim,df.attrs["ndim"]],dtype=np.float64)
    ref   = np.zeros([df.attrs["ngridmax"],twotondim],dtype=np.bool)

    iprog = 1
    istep = 10
    ncells_tot = 0

    # Loop over the cpus and read the AMR and HYDRO files in binary format
    for k in range(df.attrs["ncpu"]):

        # Print progress
        percentage = int(float(k)*100.0/float(df.attrs["ncpu"]))
        if percentage >= iprog*istep:
            print("%3i%% : read %10i cells" % (percentage,ncells_tot))
            iprog += 1

        # Read binary AMR file
        amr_fname = generate_fname(nout,path,ftype="amr",cpuid=k+1)
        with open(amr_fname, mode='rb') as amr_file:
            amrContent = amr_file.read()

        # Read binary HYDRO file
        hydro_fname = generate_fname(nout,path,ftype="hydro",cpuid=k+1)
        with open(hydro_fname, mode='rb') as hydro_file:
            hydroContent = hydro_file.read()

        # Read binary GRAVITY file
        if gravity:
            grav_fname = generate_fname(nout,path,ftype="grav",cpuid=k+1)
            with open(grav_fname, mode='rb') as grav_file:
                gravContent = grav_file.read()

        # Read binary RT file
        if rt:
            rt_fname = generate_fname(nout,path,ftype="rt",cpuid=k+1)
            with open(rt_fname, mode='rb') as rt_file:
                rtContent = rt_file.read()

        ninteg = nfloat = nlines = nstrin = nquadr = nlongi = 0

        # Need to extract info from the file header on the first loop
        if k == 0:

            # nx,ny,nz
            ninteg = 2
            nlines = 2
            [nx,ny,nz] = eng.get_binary_data(fmt="3i",content=amrContent,ninteg=ninteg,nlines=nlines)
            ncoarse = nx*ny*nz
            xbound = [float(int(nx/2)),float(int(ny/2)),float(int(nz/2))]

            # nboundary
            ninteg = 7
            nlines = 5
            [nboundary] = eng.get_binary_data(fmt="i",content=amrContent,ninteg=ninteg,nlines=nlines)
            ngridlevel = np.zeros([df.attrs["ncpu"]+nboundary,df.attrs["levelmax"]],dtype=np.int32)

            # noutput
            ninteg = 9
            nfloat = 1
            nlines = 8
            [noutput] = eng.get_binary_data(fmt="i",content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

            # dtold, dtnew
            ninteg = 12
            nfloat = 2+2*noutput
            nlines = 12
            df.attrs["dtold"] = eng.get_binary_data(fmt="%id"%(df.attrs["levelmax"]),\
                                 content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
            nfloat += 1+df.attrs["levelmax"]
            nlines += 1
            df.attrs["dtnew"] = eng.get_binary_data(fmt="%id"%(df.attrs["levelmax"]),\
                                 content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

            # hydro gamma
            ninteg = 5
            nfloat = 0
            nlines = 5
            [df.attrs["gamma"]] = eng.get_binary_data(fmt="d",content=hydroContent,ninteg=ninteg,nlines=nlines)

        # Read the number of grids
        ninteg = 14+(2*df.attrs["ncpu"]*df.attrs["levelmax"])
        nfloat = 18+(2*noutput)+(2*df.attrs["levelmax"])
        nlines = 21
        ngridlevel[:df.attrs["ncpu"],:] = np.asarray(eng.get_binary_data(fmt="%ii"%(df.attrs["ncpu"]*df.attrs["levelmax"]),\
             content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(df.attrs["levelmax"],df.attrs["ncpu"]).T

        # Read boundary grids if any
        if nboundary > 0:
            ninteg = 14+(3*df.attrs["ncpu"]*df.attrs["levelmax"])+(10*df.attrs["levelmax"])+(2*nboundary*df.attrs["levelmax"])
            nfloat = 18+(2*noutput)+(2*df.attrs["levelmax"])
            nlines = 25
            ngridlevel[df.attrs["ncpu"]:df.attrs["ncpu"]+nboundary,:] = np.asarray(eng.get_binary_data(fmt="%ii"%(nboundary*df.attrs["levelmax"]),\
                                            content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(df.attrs["levelmax"],nboundary).T

        # Determine bound key precision
        ninteg = 14+(3*df.attrs["ncpu"]*df.attrs["levelmax"])+(10*df.attrs["levelmax"])+(3*nboundary*df.attrs["levelmax"])+5
        nfloat = 18+(2*noutput)+(2*df.attrs["levelmax"])
        nlines = 21+2+3*min(1,nboundary)+1+1
        nstrin = 128
        [key_size] = eng.get_binary_data(fmt="i",content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat,nstrin=nstrin,correction=-4)

        # Offset for AMR
        ninteg1 = 14+(3*df.attrs["ncpu"]*df.attrs["levelmax"])+(10*df.attrs["levelmax"])+(3*nboundary*df.attrs["levelmax"])+5+3*ncoarse
        nfloat1 = 18+(2*noutput)+(2*df.attrs["levelmax"])
        nlines1 = 21+2+3*min(1,nboundary)+1+1+1+3
        nstrin1 = 128 + key_size

        # Offset for HYDRO
        if hydro:
            ninteg2 = 5
            nfloat2 = 1
            nlines2 = 6
            nstrin2 = 0

        # Offset for GRAVITY
        if gravity:
            ninteg3 = 4
            nfloat3 = 0
            nlines3 = 4
            nstrin3 = 0

        # Offset for RT
        if rt:
            ninteg4 = 5
            nfloat4 = 1
            nlines4 = 6
            nstrin4 = 0

        # Loop over levels
        for ilevel in range(lmax):

            # Geometry
            dxcell=0.5**(ilevel+1)
            dx2=0.5*dxcell
            for ind in range(twotondim):
                iz=int((ind)/4)
                iy=int((ind-4*iz)/2)
                ix=int((ind-2*iy-4*iz))
                xcent[ind,0]=(float(ix)-0.5)*dxcell
                xcent[ind,1]=(float(iy)-0.5)*dxcell
                xcent[ind,2]=(float(iz)-0.5)*dxcell

            # Cumulative offsets in AMR file
            ninteg_amr = ninteg1
            nfloat_amr = nfloat1
            nlines_amr = nlines1
            nstrin_amr = nstrin1

            # Cumulative offsets in HYDRO file
            if hydro:
                ninteg_hydro = ninteg2
                nfloat_hydro = nfloat2
                nlines_hydro = nlines2
                nstrin_hydro = nstrin2

            # Cumulative offsets in GRAVITY file
            if gravity:
                ninteg_grav = ninteg3
                nfloat_grav = nfloat3
                nlines_grav = nlines3
                nstrin_grav = nstrin3

            # Cumulative offsets in RT file
            if rt:
                ninteg_rt = ninteg4
                nfloat_rt = nfloat4
                nlines_rt = nlines4
                nstrin_rt = nstrin4

            # Loop over domains
            for j in range(nboundary+df.attrs["ncpu"]):

                ncache = ngridlevel[j,ilevel]

                # Skip two lines of integers
                if hydro:
                    nlines_hydro += 2
                    ninteg_hydro += 2
                if gravity:
                    nlines_grav += 2
                    ninteg_grav += 2
                if rt:
                    nlines_rt += 2
                    ninteg_rt += 2

                if ncache > 0:

                    if j == k:
                        # xg: grid coordinates
                        ninteg = ninteg_amr + ncache*3
                        nfloat = nfloat_amr
                        nlines = nlines_amr + 3
                        nstrin = nstrin_amr
                        for n in range(df.attrs["ndim"]):
                            offset = 4*ninteg + 8*(nlines+nfloat+n*(ncache+1)) + nstrin + 4
                            xg[:ncache,n] = struct.unpack("%id"%(ncache), amrContent[offset:offset+8*ncache])

                        # son indices
                        ninteg = ninteg_amr + ncache*(4+2*df.attrs["ndim"])
                        nfloat = nfloat_amr + ncache*df.attrs["ndim"]
                        nlines = nlines_amr + 4 + 3*df.attrs["ndim"]
                        nstrin = nstrin_amr
                        for ind in range(twotondim):
                            offset = 4*(ninteg+ind*ncache) + 8*(nlines+nfloat+ind) + nstrin + 4
                            son[:ncache,ind] = struct.unpack("%ii"%(ncache), amrContent[offset:offset+4*ncache])
                            # var: hydro variables
                            jvar = 0
                            if hydro:
                                for ivar in range(df.attrs["nvar_hydro"]):
                                    if var_read[ivar]:
                                        offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*df.attrs["nvar_hydro"]+ivar)*(ncache+1)) + nstrin_hydro + 4
                                        var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                        jvar += 1
                            # var: grav variables
                            if gravity:
                                for ivar in range(df.attrs["nvar_grav"]):
                                    if var_read[ivar+df.attrs["nvar_hydro"]]:
                                        offset = 4*ninteg_grav + 8*(nlines_grav+nfloat_grav+(ind*df.attrs["nvar_grav"]+ivar)*(ncache+1)) + nstrin_grav + 4
                                        var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), gravContent[offset:offset+8*ncache])
                                        jvar += 1
                            # var: rt variables
                            if rt:
                                for ivar in range(df.attrs_rt["nvar_rt"]):
                                    if var_read[ivar+df.attrs["nvar_hydro"]+df.attrs["nvar_grav"]]:
                                        offset = 4*ninteg_rt + 8*(nlines_rt+nfloat_rt+(ind*df.attrs_rt["nRTvar"]+ivar)*(ncache+1)) + nstrin_rt + 4
                                        var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), rtContent[offset:offset+8*ncache])
                                        jvar += 1
                            # var: coordinates and cell sizes
                            var[:ncache,ind,-6] = float(ilevel+1)
                            for n in range(df.attrs["ndim"]):
                                xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                var[:ncache,ind,-5+n] = xyz[:ncache,ind,n]*df.attrs["boxlen"]
                            var[:ncache,ind,-2] = dxcell*df.attrs["boxlen"]
                            var[:ncache,ind,-1] = k+1
                            # ref: True if the cell is unrefined
                            ref[:ncache,ind] = np.logical_not(np.logical_and(son[:ncache,ind] > 0, ilevel < lmax-1))

                        # Select only the cells that are in the region of interest
                        if df.attrs["ndim"] == 1:
                            cube = np.where(np.logical_and(ref[:ncache,:], \
                                            np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                                           (xyz[:ncache,:,0]-dx2)<=xmax)))
                        elif df.attrs["ndim"] == 2:
                            cube = np.where(np.logical_and(ref[:ncache,:], \
                                            np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                            np.logical_and((xyz[:ncache,:,1]+dx2)>=ymin, \
                                            np.logical_and((xyz[:ncache,:,0]-dx2)<=xmax, \
                                                           (xyz[:ncache,:,1]-dx2)<=ymax)))))
                        elif df.attrs["ndim"] == 3:
                            cube = np.where(np.logical_and(ref[:ncache,:], \
                                            np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                            np.logical_and((xyz[:ncache,:,1]+dx2)>=ymin, \
                                            np.logical_and((xyz[:ncache,:,2]+dx2)>=zmin, \
                                            np.logical_and((xyz[:ncache,:,0]-dx2)<=xmax, \
                                            np.logical_and((xyz[:ncache,:,1]-dx2)<=ymax, \
                                                           (xyz[:ncache,:,2]-dx2)<=zmax)))))))
                        else:
                            print("Bad number of dimensions")
                            return 0

                        cells = var[cube]
                        ncells = np.shape(cells)[0]
                        if ncells > 0:
                            ncells_tot += ncells
                            npieces += 1
                            # Add the cells in the master dictionary
                            data_pieces["piece"+str(npieces)] = cells

                    # Now increment the offsets while looping through the domains
                    ninteg_amr += ncache*(4+3*twotondim+2*df.attrs["ndim"])
                    nfloat_amr += ncache*df.attrs["ndim"]
                    nlines_amr += 4 + 3*twotondim + 3*df.attrs["ndim"]

                    if hydro:
                        nfloat_hydro += ncache*twotondim*df.attrs["nvar_hydro"]
                        nlines_hydro += twotondim*df.attrs["nvar_hydro"]

                    if gravity:
                        nfloat_grav += ncache*twotondim*(df.attrs["nvar_grav"])
                        nlines_grav += twotondim*(df.attrs["nvar_grav"])

                    if rt:
                        nfloat_rt += ncache*twotondim*df.attrs_rt["nvar_rt"]
                        nlines_rt += twotondim*df.attrs_rt["nvar_rt"]

            # Now increment the offsets while looping through the levels
            ninteg1 = ninteg_amr
            nfloat1 = nfloat_amr
            nlines1 = nlines_amr
            nstrin1 = nstrin_amr

            if hydro:
                ninteg2 = ninteg_hydro
                nfloat2 = nfloat_hydro
                nlines2 = nlines_hydro
                nstrin2 = nstrin_hydro

            if gravity:
                ninteg3 = ninteg_grav
                nfloat3 = nfloat_grav
                nlines3 = nlines_grav
                nstrin3 = nstrin_grav

            if rt:
                ninteg4 = ninteg_rt
                nfloat4 = nfloat_rt
                nlines4 = nlines_rt
                nstrin4 = nstrin_rt

        # # Now read particles: they are not in the loop over levels, only the cpu loop
        # if particles:
        #     fmt_to_bytes = {"b": 1 , "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8}
        #     # Read binary PARTICLE file
        #     part_fname = self.generate_fname(nout,path,ftype="part",cpuid=k+1)
        #     with open(part_fname, mode='rb') as part_file:
        #         partContent = part_file.read()
        #     # Get number of particles for this cpu
        #     offset = (fmt_to_bytes["i"] + fmt_to_bytes["e"]) * 2
        #     [npart] = eng.get_binary_data(fmt="i",content=partContent,offset=offset)
        #     if npart > 0:
        #         npart_count += npart
        #         part = np.zeros([npart, len(part_vars)],dtype=np.float64)
        #         # Determine size of localseed array
        #         offset = (fmt_to_bytes["i"] + fmt_to_bytes["e"]) * 3
        #         [recordlength] = eng.get_binary_data(fmt="i",content=partContent,offset=offset,correction=-4)
        #         localseedsize = recordlength//4
        #         # Now set offsets
        #         offset = fmt_to_bytes["i"]*(5+localseedsize) + fmt_to_bytes["e"]*8 + fmt_to_bytes["d"]*2
        #         # Go through all the particle fields and unpack the data
        #         for n in range(len(part_vars)):
        #             part[:, n] = eng.get_binary_data(fmt=("%i"%npart)+part_type[n],content=partContent,offset=offset)
        #             offset += fmt_to_bytes["e"] + npart*fmt_to_bytes[part_type[n]]

        #         # Add the cells in the master dictionary
        #         npieces_part += 1
        #         part_pieces["piece"+str(npieces_part)] = part
        # # End of reading particles ==================================================

    # Merge all the data pieces into the master data array
    master_data_array = np.concatenate(list(data_pieces.values()), axis=0)

    for i, key in enumerate(list_vars):
        df[key.strip()] = master_data_array[:, i]

    # if particles:
    #     df.attrs["npart_tot"] = npart_count
    #     master_part_array = np.concatenate(list(part_pieces.values()), axis=0)

    # Free memory
    del data_pieces,xcent,xg,son,var,xyz,ref
    # if particles:
    #     del part_pieces,part

    print("Total number of cells loaded: %i" % ncells_tot)
    # if particles:
    #     print("Total number of particles loaded: %i" % df.attrs["npart_tot"])
    # if df.attrs["nsinks"] > 0:
    #     print(("Read %i sink particle" % df.attrs["nsinks"]) + ("s" if df.attrs["nsinks"] > 1 else ""))
    # print("Generating data structure... please wait")

    # Store the number of cells
    df.attrs["ncells"] = ncells_tot

    # # We load the master array into the data structure adding one
    # # 'new_field' per variable we have read in.
    # for i in range(len(list_vars)):
    #     theKey = list_vars[i]
    #     [norm,uu] = self.get_units(theKey,df.attrs["unit_d"],df.attrs["unit_l"],df.attrs["unit_t"],df.attrs["scale"])
    #     # Replace "_" with " " to avoid error with latex when saving figures
    #     theLabel = theKey.replace("_"," ")
    #     # Use the 'new_field' function to create data field
    #     self.new_field(name=theKey,unit=uu,label=theLabel,values=(master_data_array[:,i])*norm,\
    #                    verbose=False,norm=norm,update=update,group=var_group[i])

    # # Now add new field for particles =====================================================================
    # if particles:
    #     for i in range(len(part_vars)):
    #         theKey = part_vars[i]
    #         [norm,uu] = self.get_units(theKey,df.attrs["unit_d"],df.attrs["unit_l"],df.attrs["unit_t"],df.attrs["scale"])
    #         # Replace "_" with " " to avoid error with latex when saving figures
    #         theLabel = theKey.replace("_"," ")
    #         # Use the 'new_field' function to create data field
    #         self.new_field(name=theKey,unit=uu,label=theLabel,values=master_part_array[:,i]*norm,\
    #                        verbose=False,norm=norm,update=update,group="part")

    # # Re-center the mesh around chosen center
    # self.re_center()

    return df


