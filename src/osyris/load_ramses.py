# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
import struct
import glob
from . import config as conf
from . import engine as eng

divider = "============================================"

#=======================================================================================
# This is the class which will hold the data that you read from the Ramses output
#=======================================================================================
class RamsesData(eng.OsyrisData):

    # This is the constructor which creates a `RamsesData` object.
    #
    # List of arguments and default values:
    #
    # * `nout`: (*integer*) Number of output to be read. -1 means read in the last output
    #  in the current directory. Default is `conf.default_values["nout"]`.
    #
    # * `lmax`: (*integer*) Maximum level to read up to. Default is `conf.default_values["lmax"]`.
    #
    # * `center`: (*array of 3 floats* **or** *string*) Use this to center the data coordinates
    #  around a chosen point. Possible options are a set of 3D xyz coordinates (from 0 to 1),
    #  e.g. `[0.5,0.4,0.6]`, around the density maximum `"max:density"`, around the barycentre
    #  of all cells with temperature above 1000K `"av:temperature>1000"`, or around a sink
    #  particle `"sink2"`. Default is `conf.default_values["center"]`.
    #
    # * `dx`: (*float*) Size of the region in x around the center to be read in, in units of
    #  `scale`. Default is `conf.default_values["dx"]`.
    #
    # * `dy`: (*float*) Size of the region in y around the center to be read in, in units of
    #  `scale`. Default is `conf.default_values["dy"]`.
    #
    # * `dz`: (*float*) Size of the region in z around the center to be read in, in units of
    #  `scale`. Default is `conf.default_values["dz"]`.
    #
    # * `scale`: (*float*) Spatial scale to be used for coordinates and cell sizes. Possible
    #  options are `"cm"`, `"au"` and `"pc"`. Default is `conf.default_values["scale"]`.
    #
    # * `verbose`: (*logical*) Print information about data that was read in if `True`.
    #  Default is `conf.default_values["verbose"]`.
    #
    # * `path`: (*string*) Path to the directory where to read outputs from, if different
    #  from the current directory. Default is `conf.default_values["path"]`.
    #
    # * `variables`: (*array of strings*) List of variables to be read in. To read in only
    #  the gas density and temperature, use `variables=["density","temperature"]`. Note that
    #  xyz coordinates are always read in. Default is `conf.default_values["variables"]`.
    def __init__(self,nout=conf.default_values["nout"],lmax=conf.default_values["lmax"],\
                 center=conf.default_values["center"],dx=conf.default_values["dx"],\
                 dy=conf.default_values["dy"],dz=conf.default_values["dz"],\
                 scale=conf.default_values["scale"],verbose=conf.default_values["verbose"],\
                 path=conf.default_values["path"],variables=conf.default_values["variables"]):

        # Load the Ramses data using the loader function
        status = self.data_loader(nout=nout,lmax=lmax,center=center,dx=dx,dy=dy,dz=dz,scale=scale,\
                 path=path,variables=variables)

        if status == 0:
            return

        # Convert vector components to vector containers
        self.create_vector_containers()

        # Read in custom variables if any from the configuration file
        conf.additional_variables(self)

        # Print exit message
        [var_list,typ_list] = self.get_var_list(types=True)
        print("Memory used: %.2f MB" % (typ_list.count('scalar')*self.info["ncells"]*8.0/1.0e6))
        print(self.info["infile"]+" successfully loaded")
        if verbose:
            self.print_info()
        print(divider)

        return

    #=======================================================================================

    # This function creates various file names for Ramses data.
    #
    # List of arguments and default values:
    #
    #* `nout`: (*integer*) The output number to be read in. No default value.
    #
    #* `path`: (*string*) Path to the directory where file is to be read from. Default is empty.
    #
    #* `ftype`: (*string*) The type of file. This is usually the prefix to the file, such as `amr`, `hydro` or `grav`. Default is empty.
    #
    #* `cpuid`: (*integer*) The cpu number of the file. Default is 1.
    #
    #* `ext`: (*string*) Extension for the file which is added at the end of the file name. Default is empty.
    #
    # Returns:
    #
    #* `infile`: (*string*) A file name of the form: "path/output_00001/amr_00071.out00001"+ext
    def generate_fname(self,nout,path="",ftype="",cpuid=1,ext=""):

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

    #=======================================================================================

    # This function reads in the binary Ramses output.
    #
    # List of arguments and default values:
    #
    #* `nout`: (*integer*) The output number to be read in. Default is 1.
    #
    #* `lmax`: (*integer*) Maximum level to read up to. 0 means read everything. Default is 0.
    #
    #* `center`: (*array of 3 floats* **or** *string*) Use this to center the data coordinates around a chosen point. Possible
    # options are a set of 3D xyz coordinates (from 0 to 1), e.g. `[0.5,0.4,0.6]`, around
    # the density maximum `"max:density"`, around the barycentre of all cells with temperature
    # above 1000K `"av:temperature>1000"`, or around a sink particle `"sink2"`. Default is None.
    #
    #* `dx`: (*float*) Size of the region in x around the center to be read in, in units of `scale`.
    # 0 means read everything. Default is 0.0.
    #
    #* `dy`: (*float*) Size of the region in y around the center to be read in, in units of `scale`.
    # 0 means read everything. Default is 0.0.
    #
    #* `dz`: (*float*) Size of the region in z around the center to be read in, in units of `scale`.
    # 0 means read everything. Default is 0.0.
    #
    #* `scale`: (*string*) Spatial scale to be used for coordinates and cell sizes. Possible options
    # are `"cm"`, `"au"` and `"pc"`. Default is empty.
    #
    #* `path`: (*string*) Path to the directory where to read outputs from, if different from the
    # current directory. Default is empty.
    #
    #* `variables`: (*array of strings*) List of variables to be read in. To read in only the gas density and
    # temperature, use `variables=["density","temperature"]`. Note that xyz coordinates
    # are always read in. Empty array means read everything. Default is empty.
    #
    #* `verbose`: (*logical*) Print information about data that was read in if `True` Default is False.
    #
    # Returns:
    #
    #* `status`: (*integer*) 1 is successful, 0 if not.
    def data_loader(self,nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="cm",path="",\
                    update=False,variables=[]):

        # Generate directory name from output number
        infile = self.generate_fname(nout,path)

        # Read info file and create info dictionary
        infofile = infile+"/info_"+infile.split("_")[-1]+".txt"
        status = self.read_parameter_file(fname=infofile,dict_name="info",verbose=True)
        if status < 1:
            return 0

        # Add additional information
        self.info["center"   ] = center
        self.info["scale"    ] = scale
        self.info["infile"   ] = infile
        self.info["path"     ] = path
        self.info["boxsize"  ] = self.info["boxlen"]*self.info["unit_l"]
        self.info["time"     ] = self.info["time"]*self.info["unit_t"]
        self.info["dx_load"  ] = dx
        self.info["dy_load"  ] = dy
        self.info["dz_load"  ] = dz
        self.info["lmax"     ] = lmax
        self.info["variables"] = variables
        if nout==-1:
            self.info["nout" ] = int(infile.split("_")[-1])
        else:
            self.info["nout" ] = nout

        # Read namelist file and create namelist dictionary
        nmlfile = infile+"/namelist.txt"
        self.read_parameter_file(fname=nmlfile,dict_name="namelist",evaluate=False)

        print(divider)

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
        hydro = True
        self.info["nvar_hydro"] = 0
        try:
            with open(hydrofile) as f:
                content = f.readlines()
        except IOError:
            hydro = False
        if hydro:
            # Store the total number of hydro variables
            self.info["nvar_hydro"] = len(content) - 2
            # Now add to the list of variables to be read
            for line in content[2:]:
                sp = line.split(",")
                v = sp[1].strip()
                t = sp[2].strip()
                if (len(variables) == 0) or (v in variables) or ("hydro" in variables):
                    var_read.append(True)
                    list_vars.append(v)
                    var_group.append("hydro")
                    var_type.append(t)
                else:
                    var_read.append(False)

        # Now for gravity ==================================

        # Check if self-gravity files exist
        grav_fname = self.generate_fname(nout,path,ftype="grav",cpuid=1)
        gravity = True
        self.info["nvar_grav"] = 0
        try:
            with open(grav_fname, mode='rb') as grav_file:
                gravContent = grav_file.read()
        except IOError:
            gravity = False

        # Add gravity fields
        if gravity:
            self.info["nvar_grav"] = 4
            content = ["grav_potential"]
            for n in range(self.info["ndim"]):
                content.append("grav_acceleration_"+xyz_strings[n])

            # Now add to the list of variables to be read
            for line in content:
                if (len(variables) == 0) or (line.strip() in variables) or ("gravity" in variables) or ("grav" in variables):
                    var_read.append(True)
                    list_vars.append(line.strip())
                    var_group.append("grav")
                else:
                    var_read.append(False)

        # Now for rt ==================================

        rtfile = infile+"/rt_file_descriptor.txt"
        rt = True
        self.info["nvar_rt"] = 0
        try:
            with open(rtfile) as f:
                content = f.readlines()
            # f.close()
        except IOError:
            rt = False
        if rt:
            # Store the total number of rt variables
            self.info["nvar_rt"] = len(content) - 2
            # Now add to the list of variables to be read
            for line in content[2:]:
                sp = line.split(",")
                v = sp[1].strip()
                t = sp[2].strip()
                if (len(variables) == 0) or (v in variables) or ("rt" in variables):
                    var_read.append(True)
                    list_vars.append(v)
                    var_group.append("rt")
                    var_type.append(t)
                else:
                    var_read.append(False)

        # Make sure we always read the coordinates
        var_amr = ["level","x","y","z","dx","cpu","leaf"]
        list_vars.extend(var_amr)
        var_read.extend([True] * len(var_amr))
        var_group.extend(["amr"] * len(var_amr))
        nvar_read = len(list_vars)

        # Now for particles ==================================

        # TODO: refactor this code to use the same function for reading hydro,
        # rt, and part file descriptors

        particles = True
        partfile = infile+"/part_file_descriptor.txt"
        self.info["npart_tot"] = 0
        try:
            with open(partfile) as f:
                content = f.readlines()
            # f.close()
        except IOError:
            particles = False
        if particles:
            part_read = []
            part_vars = []
            part_type = []
            # Store the total number of part variables
            self.info["nvar_part"] = len(content) - 2
            # Now add to the list of variables to be read
            for line in content[2:]:
                sp = line.split(",")
                v = "part_" + sp[1].strip()
                t = sp[2].strip()
                if (len(variables) == 0) or (v in variables) or ("part" in variables):
                    part_read.append(True)
                    part_vars.append(v)
                    part_type.append(t)
                else:
                    part_read.append(False)

        # Load sink particles if any
        self.read_sinks()

        # Find the center
        xc,yc,zc = self.find_center(dx,dy,dz)

        # Now read the amr and hydro files =============================================
        # We have to open the files in binary format, and count all the bytes in the ===
        # file structure to extract just the data we need. =============================
        # See output_amr.f90 and output_hydro.f90 in the RAMSES source. ================
        print("Processing %i files in " % (self.info["ncpu"]) + infile)

        # Define the size of the region to be read
        lconvert = conf.constants[scale]/(self.info["boxlen"]*self.info["unit_l"])
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
           lmax = self.info["levelmax"]

        # We will store the cells in a dictionary which we build as we go along.
        # The final concatenation into a single array will be done once at the end.
        data_pieces = dict()
        npieces = 0
        part_pieces = dict()
        npieces_part = 0
        npart_count = 0

        # Allocate work arrays
        twotondim = 2**self.info["ndim"]
        xcent = np.zeros([8,3],dtype=np.float64)
        xg    = np.zeros([self.info["ngridmax"],3],dtype=np.float64)
        son   = np.zeros([self.info["ngridmax"],twotondim],dtype=np.int32)
        var   = np.zeros([self.info["ngridmax"],twotondim,nvar_read],dtype=np.float64)
        xyz   = np.zeros([self.info["ngridmax"],twotondim,self.info["ndim"]],dtype=np.float64)
        ref   = np.zeros([self.info["ngridmax"],twotondim],dtype=np.bool)

        iprog = 1
        istep = 10
        ncells_tot = 0

        # Loop over the cpus and read the AMR and HYDRO files in binary format
        for k in range(self.info["ncpu"]):

            # Print progress
            percentage = int(float(k)*100.0/float(self.info["ncpu"]))
            if percentage >= iprog*istep:
                print("%3i%% : read %10i cells" % (percentage,ncells_tot))
                iprog += 1

            # Read binary AMR file
            amr_fname = self.generate_fname(nout,path,ftype="amr",cpuid=k+1)
            with open(amr_fname, mode='rb') as amr_file:
                amrContent = amr_file.read()

            # Read binary HYDRO file
            hydro_fname = self.generate_fname(nout,path,ftype="hydro",cpuid=k+1)
            with open(hydro_fname, mode='rb') as hydro_file:
                hydroContent = hydro_file.read()

            # Read binary GRAVITY file
            if gravity:
                grav_fname = self.generate_fname(nout,path,ftype="grav",cpuid=k+1)
                with open(grav_fname, mode='rb') as grav_file:
                    gravContent = grav_file.read()

            # Read binary RT file
            if rt:
                rt_fname = self.generate_fname(nout,path,ftype="rt",cpuid=k+1)
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
                ngridlevel = np.zeros([self.info["ncpu"]+nboundary,self.info["levelmax"]],dtype=np.int32)

                # noutput
                ninteg = 9
                nfloat = 1
                nlines = 8
                [noutput] = eng.get_binary_data(fmt="i",content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

                # dtold, dtnew
                ninteg = 12
                nfloat = 2+2*noutput
                nlines = 12
                self.info["dtold"] = eng.get_binary_data(fmt="%id"%(self.info["levelmax"]),\
                                     content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
                nfloat += 1+self.info["levelmax"]
                nlines += 1
                self.info["dtnew"] = eng.get_binary_data(fmt="%id"%(self.info["levelmax"]),\
                                     content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)

                # hydro gamma
                ninteg = 5
                nfloat = 0
                nlines = 5
                [self.info["gamma"]] = eng.get_binary_data(fmt="d",content=hydroContent,ninteg=ninteg,nlines=nlines)

            # Read the number of grids
            ninteg = 14+(2*self.info["ncpu"]*self.info["levelmax"])
            nfloat = 18+(2*noutput)+(2*self.info["levelmax"])
            nlines = 21
            ngridlevel[:self.info["ncpu"],:] = np.asarray(eng.get_binary_data(fmt="%ii"%(self.info["ncpu"]*self.info["levelmax"]),\
                 content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(self.info["levelmax"],self.info["ncpu"]).T

            # Read boundary grids if any
            if nboundary > 0:
                ninteg = 14+(3*self.info["ncpu"]*self.info["levelmax"])+(10*self.info["levelmax"])+(2*nboundary*self.info["levelmax"])
                nfloat = 18+(2*noutput)+(2*self.info["levelmax"])
                nlines = 25
                ngridlevel[self.info["ncpu"]:self.info["ncpu"]+nboundary,:] = np.asarray(eng.get_binary_data(fmt="%ii"%(nboundary*self.info["levelmax"]),\
                                                content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)).reshape(self.info["levelmax"],nboundary).T

            # Determine bound key precision
            ninteg = 14+(3*self.info["ncpu"]*self.info["levelmax"])+(10*self.info["levelmax"])+(3*nboundary*self.info["levelmax"])+5
            nfloat = 18+(2*noutput)+(2*self.info["levelmax"])
            nlines = 21+2+3*min(1,nboundary)+1+1
            nstrin = 128
            [key_size] = eng.get_binary_data(fmt="i",content=amrContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat,nstrin=nstrin,correction=-4)

            # Offset for AMR
            ninteg1 = 14+(3*self.info["ncpu"]*self.info["levelmax"])+(10*self.info["levelmax"])+(3*nboundary*self.info["levelmax"])+5+3*ncoarse
            nfloat1 = 18+(2*noutput)+(2*self.info["levelmax"])
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
                for j in range(nboundary+self.info["ncpu"]):

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
                            for n in range(self.info["ndim"]):
                                offset = 4*ninteg + 8*(nlines+nfloat+n*(ncache+1)) + nstrin + 4
                                xg[:ncache,n] = struct.unpack("%id"%(ncache), amrContent[offset:offset+8*ncache])

                            # son indices
                            ninteg = ninteg_amr + ncache*(4+2*self.info["ndim"])
                            nfloat = nfloat_amr + ncache*self.info["ndim"]
                            nlines = nlines_amr + 4 + 3*self.info["ndim"]
                            nstrin = nstrin_amr
                            for ind in range(twotondim):
                                offset = 4*(ninteg+ind*ncache) + 8*(nlines+nfloat+ind) + nstrin + 4
                                son[:ncache,ind] = struct.unpack("%ii"%(ncache), amrContent[offset:offset+4*ncache])
                                # var: hydro variables
                                jvar = 0
                                if hydro:
                                    for ivar in range(self.info["nvar_hydro"]):
                                        if var_read[ivar]:
                                            offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*self.info["nvar_hydro"]+ivar)*(ncache+1)) + nstrin_hydro + 4
                                            var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                            jvar += 1
                                # var: grav variables
                                if gravity:
                                    for ivar in range(self.info["nvar_grav"]):
                                        if var_read[ivar+self.info["nvar_hydro"]]:
                                            offset = 4*ninteg_grav + 8*(nlines_grav+nfloat_grav+(ind*self.info["nvar_grav"]+ivar)*(ncache+1)) + nstrin_grav + 4
                                            var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), gravContent[offset:offset+8*ncache])
                                            jvar += 1
                                # var: rt variables
                                if rt:
                                    for ivar in range(self.info_rt["nvar_rt"]):
                                        if var_read[ivar+self.info["nvar_hydro"]+self.info["nvar_grav"]]:
                                            offset = 4*ninteg_rt + 8*(nlines_rt+nfloat_rt+(ind*self.info_rt["nRTvar"]+ivar)*(ncache+1)) + nstrin_rt + 4
                                            var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), rtContent[offset:offset+8*ncache])
                                            jvar += 1
                                # var: coordinates and cell sizes
                                var[:ncache,ind,-7] = float(ilevel+1)
                                for n in range(self.info["ndim"]):
                                    xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                    var[:ncache,ind,-6+n] = xyz[:ncache,ind,n]*self.info["boxlen"]
                                var[:ncache,ind,-3] = dxcell*self.info["boxlen"]
                                var[:ncache,ind,-2] = k+1
                                # leaf cells: True if the cell is unrefined
                                var[:ncache,ind,-1] = 1.0 * np.logical_not(np.logical_and(son[:ncache,ind] > 0,ilevel < lmax-1))

                            # Select only the cells that are in the region of interest
                            if self.info["ndim"] == 1:
                                cube = np.where(np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                                               (xyz[:ncache,:,0]-dx2)<=xmax))
                            elif self.info["ndim"] == 2:
                                cube = np.where(np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                                np.logical_and((xyz[:ncache,:,1]+dx2)>=ymin, \
                                                np.logical_and((xyz[:ncache,:,0]-dx2)<=xmax, \
                                                               (xyz[:ncache,:,1]-dx2)<=ymax))))
                            elif self.info["ndim"] == 3:
                                cube = np.where(np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                                np.logical_and((xyz[:ncache,:,1]+dx2)>=ymin, \
                                                np.logical_and((xyz[:ncache,:,2]+dx2)>=zmin, \
                                                np.logical_and((xyz[:ncache,:,0]-dx2)<=xmax, \
                                                np.logical_and((xyz[:ncache,:,1]-dx2)<=ymax, \
                                                               (xyz[:ncache,:,2]-dx2)<=zmax))))))
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
                        ninteg_amr += ncache*(4+3*twotondim+2*self.info["ndim"])
                        nfloat_amr += ncache*self.info["ndim"]
                        nlines_amr += 4 + 3*twotondim + 3*self.info["ndim"]

                        if hydro:
                            nfloat_hydro += ncache*twotondim*self.info["nvar_hydro"]
                            nlines_hydro += twotondim*self.info["nvar_hydro"]

                        if gravity:
                            nfloat_grav += ncache*twotondim*(self.info["nvar_grav"])
                            nlines_grav += twotondim*(self.info["nvar_grav"])

                        if rt:
                            nfloat_rt += ncache*twotondim*self.info_rt["nvar_rt"]
                            nlines_rt += twotondim*self.info_rt["nvar_rt"]

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

            # Now read particles: they are not in the loop over levels, only the cpu loop
            if particles:
                fmt_to_bytes = {"b": 1 , "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8}
                # Read binary PARTICLE file
                part_fname = self.generate_fname(nout,path,ftype="part",cpuid=k+1)
                with open(part_fname, mode='rb') as part_file:
                    partContent = part_file.read()
                # Get number of particles for this cpu
                offset = (fmt_to_bytes["i"] + fmt_to_bytes["e"]) * 2
                [npart] = eng.get_binary_data(fmt="i",content=partContent,offset=offset)
                if npart > 0:
                    npart_count += npart
                    part = np.zeros([npart, len(part_vars)],dtype=np.float64)
                    # Determine size of localseed array
                    offset = (fmt_to_bytes["i"] + fmt_to_bytes["e"]) * 3
                    [recordlength] = eng.get_binary_data(fmt="i",content=partContent,offset=offset,correction=-4)
                    localseedsize = recordlength//4
                    # Now set offsets
                    offset = fmt_to_bytes["i"]*(5+localseedsize) + fmt_to_bytes["e"]*8 + fmt_to_bytes["d"]*2
                    # Go through all the particle fields and unpack the data
                    for n in range(len(part_vars)):
                        part[:, n] = eng.get_binary_data(fmt=("%i"%npart)+part_type[n],content=partContent,offset=offset)
                        offset += fmt_to_bytes["e"] + npart*fmt_to_bytes[part_type[n]]

                    # Add the cells in the master dictionary
                    npieces_part += 1
                    part_pieces["piece"+str(npieces_part)] = part
            # End of reading particles ==================================================

        # Merge all the data pieces into the master data array
        master_data_array = np.concatenate(list(data_pieces.values()), axis=0)
        if particles:
            self.info["npart_tot"] = npart_count
            master_part_array = np.concatenate(list(part_pieces.values()), axis=0)

        # Free memory
        del data_pieces,xcent,xg,son,var,xyz,ref
        if particles:
            del part_pieces,part

        print("Total number of cells loaded: %i" % ncells_tot)
        if particles:
            print("Total number of particles loaded: %i" % self.info["npart_tot"])
        if self.info["nsinks"] > 0:
            print(("Read %i sink particle" % self.info["nsinks"]) + ("s" if self.info["nsinks"] > 1 else ""))
        print("Generating data structure... please wait")

        # Store the number of cells
        self.info["ncells"] = ncells_tot

        # Finally we add one 'new_field' per variable we have read in =========================================
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            [norm,uu] = self.get_units(theKey,self.info["unit_d"],self.info["unit_l"],self.info["unit_t"],self.info["scale"])
            # Replace "_" with " " to avoid error with latex when saving figures
            theLabel = theKey.replace("_"," ")
            # Use the 'new_field' function to create data field
            self.new_field(name=theKey,unit=uu,label=theLabel,values=master_data_array[:,i]*norm,\
                           verbose=False,norm=norm,update=update,group=var_group[i])

        # Now add new field for particles =====================================================================
        if particles:
            for i in range(len(part_vars)):
                theKey = part_vars[i]
                [norm,uu] = self.get_units(theKey,self.info["unit_d"],self.info["unit_l"],self.info["unit_t"],self.info["scale"])
                # Replace "_" with " " to avoid error with latex when saving figures
                theLabel = theKey.replace("_"," ")
                # Use the 'new_field' function to create data field
                self.new_field(name=theKey,unit=uu,label=theLabel,values=master_part_array[:,i]*norm,\
                               verbose=False,norm=norm,update=update,group="part")

        # Finally, add some useful information to save compute time later
        self.info["levelmax_active"] = np.nanmax(self.level.values)
        self.info["leafs"] = np.where(self.leaf.values == 1.0)

        # Re-center the mesh around chosen center
        self.re_center()

        return 1

    #=======================================================================================

    # Read in sink particle `.csv` file if present.
    def read_sinks(self):

        sinkfile = self.info["infile"]+"/sink_"+self.info["infile"].split("_")[-1]+".csv"
        try:
            with open(sinkfile) as f:
                content = f.readlines()
        except IOError:
            self.info["nsinks"] = 0
            return
        # Read the file header to get information on fields
        sink_vars = content[0].rstrip().replace(" # ", "").split(",")
        sink_units = content[1].rstrip().replace(" # ", "").split(",")
        self.info["nsinks"] = len(content) - 2
        if self.info["nsinks"] > 0:
            self.sinks = dict()
            for entry in sink_vars:
                self.sinks[entry] = np.zeros(self.info["nsinks"], dtype=np.float64)
            for i in range(self.info["nsinks"]):
                line = np.asarray(content[i+2].rstrip().split(","), dtype=np.float64)
                for j, entry in enumerate(sink_vars):
                    self.sinks[entry][i] = np.float64(line[j])
            self.sinks["x"] *= self.info["unit_l"]
            self.sinks["y"] *= self.info["unit_l"]
            self.sinks["z"] *= self.info["unit_l"]
            self.sinks["id"] = np.int32(self.sinks["id"])
            self.sinks["level"] = np.int32(self.sinks["level"])
            self.sinks["radius"] = 4.0/(2.0**self.sinks["level"])

        return

    #=======================================================================================

    # This function updates all the fields of a RamsesData container with values from a new
    # output number, including derived fields.
    #
    # List of arguments and default values:
    #
    #* `nout`: (*integer*) The output number to be read in. Default is 1.
    #
    #* `lmax`: (*integer*) Maximum level to read up to. 0 means read everything. Default is 0.
    #
    #* `center`: (*array of 3 floats* **or** *string*) Use this to center the data coordinates around a chosen point.
    # Possible options are a set of 3D xyz coordinates (from 0 to 1), e.g. `[0.5,0.4,0.6]`,
    # around the density maximum `"max:density"`, around the barycentre of all cells with
    # temperature above 1000K `"av:temperature>1000"`, or around a sink particle `"sink2"`.
    # Default is None.
    #
    #* `dx`: (*float*) Size of the region in x around the center to be read in, in units of `scale`.
    # 0 means read everything. Default is 0.0.
    #
    #* `dy`: (*float*) Size of the region in y around the center to be read in, in units of `scale`.
    # 0 means read everything. Default is 0.0.
    #
    #* `dz`: (*float*) Size of the region in z around the center to be read in, in units of `scale`.
    # 0 means read everything. Default is 0.0.
    #
    #* `scale`: (*string*) Spatial scale to be used for coordinates and cell sizes. Possible options
    # are `"cm"`, `"au"` and `"pc"`. Default is empty.
    #
    #* `path`: (*string*) Path to the directory where to read outputs from, if different from the
    # current directory. Default is empty.
    #
    #* `variables`: (*array of strings*) List of variables to be read in. To read in only the gas density and
    # temperature, use `variables=["density","temperature"]`. Note that xyz coordinates
    # are always read in. Empty array means read everything. Default is empty.
    #
    #* `verbose`: (*logical*) Print information about data that was read in if `True` Default is False.
    def update_values(self,nout=-1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="",\
                      path="",variables=[],verbose=False):

        ## Check if new output number is requested. If not, use same nout as before
        #if nout == "none":
            #nout = self.info["nout"]

        # Check if new lmax is requested. If not, use same lmax as before
        if lmax == 0:
            lmax = self.info["lmax"]

        # Check if a new center is requested. If not, use same center as before
        try:
            dummy = len(center)
        except TypeError:
            center = self.info["center"]

        # Check if new scale is requested. If not, use same scale as before
        if len(scale) == 0:
            scale = self.info["scale"]

        # Check if new path is requested. If not, use same path as before
        if len(path) == 0:
            path = self.info["path"]

        # Check if new dx,dy,dz are requested. If not, use same as before
        if dx == 0.0:
            dx = self.info["dx_load"]
        if dy == 0.0:
            dy = self.info["dy_load"]
        if dz == 0.0:
            dz = self.info["dz_load"]

        # Check if new list of variables is requested. If not, use same list as before
        if len(variables) == 0:
            variables = self.info["variables"]

        # Load the Ramses data using the loader function
        status = self.data_loader(nout=nout,lmax=lmax,center=center,dx=dx,dy=dy,dz=dz,scale=scale,\
                                  path=path,variables=variables,update=True)

        if status == 0:
            return

        # Now go through the fields and update the values of fields that have an operation
        # attached to them. IMPORTANT!: this needs to be done in the right order: use the
        # depth key to determine which variables depend on others
        key_list = self.get_var_list()
        key_list = sorted(key_list,key=lambda x:getattr(self,x).depth)
        with np.errstate(divide="ignore"):
            for key in key_list:
                dataField = getattr(self,key)
                if len(dataField.operation) > 0:
                    print("Re-computing "+key)
                    dataField.values = eval(dataField.operation)

        ## Re-center the mesh around chosen center
        #self.re_center()

        print("Data successfully updated with values from "+self.info["infile"])
        if verbose:
            self.print_info()
        print(divider)

        return

    #=======================================================================================

    # This function returns the appropriate scaling for a variable which was read
    # in code units by the data loader. It tries to identify if we are dealing with a
    # density or a pressure and returns the appropriate combination of ud, ul and ut. It
    # also returns the unit as a string for plotting on the axes.
    #
    # List of arguments and default values:
    #
    #* `string`: (*string*) Name of the variable. There is no default.
    #
    #* `ud`: (*float*) Denstiy scaling. There is no default.
    #
    #* `ul`: (*float*) Length scaling. There is no default.
    #
    #* `ut`: (*float*) Time scaling. There is no default.
    #
    #* `scale`: (*string*) String to be used for scaling units on axes. Default is `"cm"`.
    #
    # Returns:
    #
    #* `norm`: (*float*) A scaling factor to convert from code units to cgs.
    #
    #* `label`: (*string*) A string describing the units, to be used on axes.
    def get_units(self,string,ud,ul,ut,scale="cm"):
        if string == "density":
            return [ud,"g/cm3"]
        elif (string.startswith("velocity")) or (string.startswith("part_velocity")):
            return [ul/ut,"cm/s"]
        elif string.startswith("momentum"):
            return [ud*ul/ut,"g/cm2/s"]
        elif (string.startswith("B_")) or (string.startswith("part_tracer_b")):
            return [np.sqrt(4.0*np.pi*ud*(ul/ut)**2),"G"]
        elif (string.startswith("current_")):
            return [np.sqrt(4.0*np.pi*ud*(ul/ut)**2)/(self.info["boxlen"]*ul),"G/cm"]
        elif ("acceleration" in string):
            return [ul/ut**2,"cm/s2"]
        elif string == ("thermal_pressure") or (string.count("energy") > 0):
            return [ud*((ul/ut)**2),"erg/cm3"]
        elif (string == "x") or (string == "y") or (string == "z") or (string == "dx"):
            return [ul,scale]
        elif string.startswith("part_position"):
            return [ul*self.info["boxlen"],scale]
        elif string == "temperature":
            return [1.0,"K"]
        elif string.startswith("photon_density"):
            return [self.info_rt["unit_np"],"photons/cm3"]
        elif string.startswith("photon_flux"):
            return [self.info_rt["unit_pf"],"photons/cm2/s"]
        else:
            for key in conf.default_units.keys():
                if string == key:
                    new_string = conf.default_units[string][0].replace("unit_d","self.info[\"unit_d\"]")
                    new_string = new_string.replace("unit_l","self.info[\"unit_l\"]")
                    new_string = new_string.replace("unit_t","self.info[\"unit_t\"]")
                    uu = eval(new_string)
                    return [uu,conf.default_units[string][1]]
            return [1.0,""]
