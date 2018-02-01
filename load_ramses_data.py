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

import numpy as np
import struct
import glob
import config_osiris as conf

divider = "============================================"



class OsirisData():
    
    def __init__(self,values=None,unit=None,label=None,operation=None,depth=None,norm=1.0,\
                 parent=None,kind='scalar',vec_x=False,vec_y=False,vec_z=False,name="",vector_component=False):
        
        self.values = values
        self.unit = unit
        self.label = label
        self.operation = operation
        self.depth = depth
        self.norm = norm
        self.kind = kind
        self.parent = parent
        self.name = name
        if vec_x:
            self.x = vec_x
        if vec_y:
            self.y = vec_y
        if vec_z:
            self.z = vec_z
        self.vector_component = vector_component
        
        return
    


#=======================================================================================
# This is the class which will hold the data that you read from the Ramses output
# It calls "rd.ramses_data" which is a fortran file reader.
# It then stores the data in a dictionary named "data"
#=======================================================================================
class LoadRamsesData():
 
    #===================================================================================
    # The constructor reads in the data and fills the data structure which is a python
    # dictionary. The arguments are:
    # - nout  : the number of the output. It can be -1 for the last output
    # - lmax  : maximum AMR level to be read
    # - center: used to re-centre the mesh coordinates around a given center. Possible
    #           values are and array of 3 numbers between 0 and 1, e.g. [0.51,0.46,0.33]
    #           or you can use center="auto" to automatically find the densest cell
    # - dx    : size of the domain to be read in the x dimension, in units of scale
    # - dy    : size of the domain to be read in the y dimension, in units of scale
    # - dz    : size of the domain to be read in the z dimension, in units of scale
    # - scale : spatial scale conversion for distances. Possible values are "cm", "au"
    #           and "pc"
    #===================================================================================
    def __init__(self,nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale=False,verbose=False,\
                 path="",variables=[]):
                
        # Load the Ramses data using the loader function
        status = self.data_loader(nout=nout,lmax=lmax,center=center,dx=dx,dy=dy,dz=dz,scale=scale,\
                 path=path,variables=variables)
        
        if status == 0:
            return
        
        ## Re-center the mesh around chosen center
        #self.re_center()
        
        # Read in custom variables if any from the configuration file
        conf.additional_variables(self)
        
        # Convert vector components to vector containers
        self.create_vector_containers()
        
        # Print exit message
        [var_list,typ_list] = self.get_var_list(types=True)
        print("Memory used: %.2f Mb" % (typ_list.count('scalar')*self.info["ncells"]*8.0/1.0e6))
        print(self.info["infile"]+" successfully loaded")
        if verbose:
            self.print_info()
        print(divider)
    
    #=======================================================================================
    # Generate the file name
    #=======================================================================================
    def generate_fname(self,nout,path="",ftype="",cpuid=1):
        
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
            infile = infile+"/"+ftype+"_"+number+".out"+str(cpuid).zfill(5)
            
        return infile
    
    #=======================================================================================
    # Load the data from fortran routine
    #=======================================================================================
    def data_loader(self,nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="cm",path="",\
                    update=False,variables=[]):
        
        # Generate filename from output number
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
        
        # Determine whether self-gravity was used and is to be read
        grav_fname = self.generate_fname(nout,path,ftype="grav",cpuid=1)
        try:
            with open(grav_fname, mode='rb') as grav_file: # b is important -> binary
                gravContent = grav_file.read()
            grav_file.close()
            gravity = True
        except IOError:
            gravity = False
        
        # Read the number of variables from the hydro_file_descriptor.txt
        # and select the ones to be read if specified by user
        hydrofile = infile+"/hydro_file_descriptor.txt"
        try:
            with open(hydrofile) as f:
                content = f.readlines()
            f.close()
        except IOError:
            # If hydro_file_descriptor.txt does not exist, mimic the
            # content by using the default names from the config file
            content = ["nvar = "+str(len(conf.default_values["var_names"]))]
            ivar = 0
            for var in conf.default_values["var_names"]:
                ivar = ivar + 1
                content.append("variable #"+str(ivar)+": "+var)
        # Read the total number of hydro variables
        for line in content:
            sp = line.split("=")
            if len(sp) > 1:
                if sp[0].strip() == "nvar":
                    self.info["nvar_hydro"] = int(sp[1].strip())
                    break
        
        # Count variables
        nv_count = self.info["nvar_hydro"]+6
        
        # Add gravity fields
        if gravity:
            xyz_strings = "xyz"
            ivar = self.info["nvar_hydro"]+1
            content.append("variable #"+str(ivar)+": grav_potential")
            for n in range(self.info["ndim"]):
                ivar += 1
                content.append("variable #"+str(ivar)+": grav_acceleration_"+xyz_strings[n])
            nv_count += 1+self.info["ndim"]
            
        # Now go through all the variables and check if they are to be read or skipped        
        var_read = np.ones([nv_count],dtype=np.bool)
        list_vars = []
        ivar = 0
        for line in content:
            sp = line.split(":")
            if len(sp) > 1:
                if (len(variables) == 0) or (variables.count(sp[1].strip()) > 0):
                    #var_read += "1 "
                    var_read[ivar] = True
                    list_vars.append(sp[1].strip())
                else:
                    var_read[ivar] = False
                ivar += 1
                
        # Make sure we always read the coordinates
        list_vars.extend(("level","x","y","z","dx","cpu"))
        nvar_read = len(list_vars)
        
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
            with open(amr_fname, mode='rb') as amr_file: # b is important -> binary
                amrContent = amr_file.read()
            amr_file.close()
            
            # Read binary HYDRO file
            hydro_fname = self.generate_fname(nout,path,ftype="hydro",cpuid=k+1)
            with open(hydro_fname, mode='rb') as hydro_file: # b is important -> binary
                hydroContent = hydro_file.read()
            hydro_file.close()
            
            # Read binary GRAVITY file
            if gravity:
                grav_fname = self.generate_fname(nout,path,ftype="grav",cpuid=k+1)
                with open(grav_fname, mode='rb') as grav_file: # b is important -> binary
                    gravContent = grav_file.read()
                grav_file.close()
            
            # Need to extract info from the file header on the first loop
            if k == 0:
            
                # nx,ny,nz
                ninteg = 2
                nfloat = 0
                nlines = 2
                nstrin = 0
                nquadr = 0
                offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
                [nx,ny,nz] = struct.unpack("3i", amrContent[offset:offset+12])
                ncoarse = nx*ny*nz
                xbound = [float(int(nx/2)),float(int(ny/2)),float(int(nz/2))]
                
                # nboundary
                ninteg = 7
                nfloat = 0
                nlines = 5
                nstrin = 0
                nquadr = 0
                offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
                nboundary = struct.unpack("i", amrContent[offset:offset+4])[0]
                ngridlevel = np.zeros([self.info["ncpu"]+nboundary,self.info["levelmax"]],dtype=np.int32)
                
                # noutput
                ninteg = 9
                nfloat = 1
                nlines = 8
                nstrin = 0
                nquadr = 0
                offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
                noutput = struct.unpack("i", amrContent[offset:offset+4])[0]
                
                # hydro gamma
                ninteg = 5
                nfloat = 0
                nlines = 5
                nstrin = 0
                offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
                self.info["gamma"] = struct.unpack("d", hydroContent[offset:offset+8])[0]
                
                # dtold, dtnew
                ninteg = 12
                nfloat = 2+2*noutput
                nlines = 12
                nstrin = 0
                nquadr = 0
                offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
                self.info["dtold"] = struct.unpack("%id"%(self.info["levelmax"]), amrContent[offset:offset+8*self.info["levelmax"]])
                nfloat = 3+2*noutput+self.info["levelmax"]
                nlines = 13
                self.info["dtnew"] = struct.unpack("%id"%(self.info["levelmax"]), amrContent[offset:offset+8*self.info["levelmax"]])

            # Read the number of grids
            ninteg = 14+(2*self.info["ncpu"]*self.info["levelmax"])
            nfloat = 18+(2*noutput)+(2*self.info["levelmax"])
            nlines = 21
            nstrin = 0
            nquadr = 0
            offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
            ngridlevel[:self.info["ncpu"],:] = np.asarray(struct.unpack("%ii"%(self.info["ncpu"]*self.info["levelmax"]), amrContent[offset:offset+4*self.info["ncpu"]*self.info["levelmax"]])).reshape(self.info["levelmax"],self.info["ncpu"]).T
            
            # Read boundary grids if any
            if nboundary > 0:
                ninteg = 14+(3*self.info["ncpu"]*self.info["levelmax"])+(10*self.info["levelmax"])+(2*nboundary*self.info["levelmax"])
                nfloat = 18+(2*noutput)+(2*self.info["levelmax"])
                nlines = 25
                nstrin = 0
                nquadr = 0
                offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
                ngridlevel[self.info["ncpu"]:self.info["ncpu"]+nboundary,:] = np.asarray(struct.unpack("%ii"%(nboundary*self.info["levelmax"]), amrContent[offset:offset+4*nboundary*self.info["levelmax"]])).reshape(self.info["levelmax"],nboundary).T
    
            # Determine bound key precision
            ninteg = 14+(3*self.info["ncpu"]*self.info["levelmax"])+(10*self.info["levelmax"])+(3*nboundary*self.info["levelmax"])+5
            nfloat = 18+(2*noutput)+(2*self.info["levelmax"])
            nlines = 21+2+3*min(1,nboundary)+1+1
            nstrin = 128
            nquadr = 0
            offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16
            key_size = struct.unpack("i", amrContent[offset:offset+4])[0]
            
            # Offset for AMR
            ninteg1 = 14+(3*self.info["ncpu"]*self.info["levelmax"])+(10*self.info["levelmax"])+(3*nboundary*self.info["levelmax"])+5+3*ncoarse
            nfloat1 = 18+(2*noutput)+(2*self.info["levelmax"])
            nlines1 = 21+2+3*min(1,nboundary)+1+1+1+3
            nstrin1 = 128 + key_size
            
            # Offset for HYDRO
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
                                
                # Loop over domains
                for j in range(nboundary+self.info["ncpu"]):
                    
                    ncache = ngridlevel[j,ilevel]
                    
                    # Skip two lines of integers
                    nlines_hydro += 2
                    ninteg_hydro += 2
                    if gravity:
                        nlines_grav += 2
                        ninteg_grav += 2
                    
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
                                for ivar in range(self.info["nvar_hydro"]):
                                    if var_read[ivar]:
                                        offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*self.info["nvar_hydro"]+ivar)*(ncache+1)) + nstrin_hydro + 4
                                        var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                        jvar += 1
                                # var: grav variables
                                if gravity:
                                    #jvar = 0
                                    for ivar in range(self.info["ndim"]+1):
                                        if var_read[ivar+self.info["nvar_hydro"]]:
                                            offset = 4*ninteg_grav + 8*(nlines_grav+nfloat_grav+(ind*(self.info["ndim"]+1)+ivar)*(ncache+1)) + nstrin_grav + 4
                                            var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), gravContent[offset:offset+8*ncache])
                                            jvar += 1
                                var[:ncache,ind,-6] = float(ilevel+1)
                                for n in range(self.info["ndim"]):
                                    xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                    var[:ncache,ind,-5+n] = xyz[:ncache,ind,n]*self.info["boxlen"]
                                var[:ncache,ind,-2] = dxcell*self.info["boxlen"]
                                var[:ncache,ind,-1] = k+1
                                # ref: True if the cell is unrefined
                                ref[:ncache,ind] = np.logical_not(np.logical_and(son[:ncache,ind] > 0,ilevel < lmax-1))

                            # Select only the unrefined cells that are in the region of interest
                            if self.info["ndim"] == 1:
                                cube = np.where(np.logical_and(ref[:ncache,:], \
                                                np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                                               (xyz[:ncache,:,0]-dx2)<=xmax)))
                            elif self.info["ndim"] == 2:
                                cube = np.where(np.logical_and(ref[:ncache,:], \
                                                np.logical_and((xyz[:ncache,:,0]+dx2)>=xmin, \
                                                np.logical_and((xyz[:ncache,:,1]+dx2)>=ymin, \
                                                np.logical_and((xyz[:ncache,:,0]-dx2)<=xmax, \
                                                               (xyz[:ncache,:,1]-dx2)<=ymax)))))
                            elif self.info["ndim"] == 3:
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
                        ninteg_amr += ncache*(4+3*twotondim+2*self.info["ndim"])
                        nfloat_amr += ncache*self.info["ndim"]
                        nlines_amr += 4 + 3*twotondim + 3*self.info["ndim"]
                        
                        nfloat_hydro += ncache*twotondim*self.info["nvar_hydro"]
                        nlines_hydro += twotondim*self.info["nvar_hydro"]
                        
                        if gravity:
                            nfloat_grav += ncache*twotondim*(self.info["ndim"]+1)
                            nlines_grav += twotondim*(self.info["ndim"]+1)
                
                # Now increment the offsets while looping through the levels
                ninteg1 = ninteg_amr
                nfloat1 = nfloat_amr
                nlines1 = nlines_amr
                nstrin1 = nstrin_amr
                
                ninteg2 = ninteg_hydro
                nfloat2 = nfloat_hydro
                nlines2 = nlines_hydro
                nstrin2 = nstrin_hydro
                
                ninteg3 = ninteg_grav
                nfloat3 = nfloat_grav
                nlines3 = nlines_grav
                nstrin3 = nstrin_grav
        
        # Merge all the data pieces into the master data array
        master_data_array = np.concatenate(list(data_pieces.values()), axis=0)
        
        # Free memory
        del data_pieces,xcent,xg,son,var,xyz,ref
        
        print("Total number of cells loaded: %i" % ncells_tot)
        if self.info["nsinks"] > 0:
            print("Read %i sink particles" % self.info["nsinks"])
        print("Generating data structure... please wait")
        
        # Store the number of cells
        self.info["ncells"] = ncells_tot
        
        # This is the master data dictionary. For each entry, the dict has 5 fields.
        # It loops through the list of variables that it got from the file loader.
        #if not update:
            #self.data = dict()
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            #if not update:
                #self.data[theKey] = dict()
            [norm,uu] = self.get_units(theKey,self.info["unit_d"],self.info["unit_l"],self.info["unit_t"],self.info["scale"])
            # Replace "_" with " " to avoid error with latex when saving figures
            theLabel = theKey.replace("_"," ")
            # Use the 'new_field' function to create data field
            self.new_field(name=theKey,operation="",unit=uu,label=theLabel,values=master_data_array[:,i]*norm,verbose=False,norm=norm,update=update)
        
                        
                        
        #self.print_info()
        
        # Hard coded additional data fields needed
        [norm,uu] = self.get_units("x",self.info["unit_d"],self.info["unit_l"],self.info["unit_t"],self.info["scale"])
        self.new_field(name="x_raw",operation="x",unit=uu,label="x raw",verbose=False,norm=norm,update=update)
        self.new_field(name="y_raw",operation="y",unit=uu,label="y raw",verbose=False,norm=norm,update=update)
        self.new_field(name="z_raw",operation="z",unit=uu,label="z raw",verbose=False,norm=norm,update=update)
        self.new_field(name="dx_raw",operation="dx",unit=uu,label="dx raw",verbose=False,norm=norm,update=update)
        self.new_field(name="x_box",values=self.get("x")/norm/self.info["boxlen"],unit="",label="x box",verbose=False,norm=1.0,update=update)
        self.new_field(name="y_box",values=self.get("y")/norm/self.info["boxlen"],unit="",label="y box",verbose=False,norm=1.0,update=update)
        self.new_field(name="z_box",values=self.get("z")/norm/self.info["boxlen"],unit="",label="z box",verbose=False,norm=1.0,update=update)
        self.new_field(name="dx_box",values=self.get("dx")/norm/self.info["boxlen"],unit="",label="dx box",verbose=False,norm=1.0,update=update)

        #self.print_info()
        
        # Re-center the mesh around chosen center
        self.re_center()

        return 1
    
    #=======================================================================================
    # Print information about the data that was loaded.
    #=======================================================================================
    def read_parameter_file(self,fname="",dict_name="",evaluate=True,verbose=False):
    
        # Read info file and create dictionary
        try:
            with open(fname) as f:
                content = f.readlines()
            f.close()
        except IOError:
            # Clean exit if the file was not found
            if verbose:
                print("File not found: "+fname)
            #if raise_error:
                #raise IOError
            #else:
            return 0
        
        setattr(self,dict_name,dict())
        for line in content:
            sp = line.split("=")
            if len(sp) > 1:
                if evaluate:
                    try:
                        getattr(self,dict_name)[sp[0].strip()] = eval(sp[1].strip())
                    except (NameError,SyntaxError):
                        getattr(self,dict_name)[sp[0].strip()] = sp[1].strip()
                else:
                    getattr(self,dict_name)[sp[0].strip()] = sp[1].strip()

        return 1
    
    #=======================================================================================
    # Print information about the data that was loaded.
    #=======================================================================================
    def print_info(self):
        print("--------------------------------------------")
        for key in sorted(self.info.keys()):
            theShape = np.shape(self.info[key])
            if len(theShape) > 0:
                try:
                    print(key+": ["+str(self.info[key][0])+" ... "+str(self.info[key][-1])+"]")
                except IndexError:
                    print(key+": "+str(self.info[key]))
            else:
                print(key+": "+str(self.info[key]))
        print("--------------------------------------------")
        maxlen1 = maxlen2 = maxlen3 = maxlen4 = maxlen5 = 0
        key_list = self.get_var_list()
        print_list = dict()
        #key_list = sorted(key_list,key=lambda x:len(x),reverse=True)
        for key in sorted(key_list):
            print_list[key] = []
            print_list[key].append(key)
            maxlen1 = max(maxlen1,len(key))
            print_list[key].append(getattr(self,key).kind)
            maxlen2 = max(maxlen2,len(print_list[key][1]))
            print_list[key].append(getattr(self,key).unit)
            maxlen3 = max(maxlen3,len(print_list[key][2]))
            if print_list[key][1] == 'vector':
                print_list[key].append("--")
                print_list[key].append("--")
            else:
                try:
                    print_list[key].append(str(np.nanmin(getattr(self,key).values)))
                except TypeError:
                    print_list[key].append("--")
                try:
                    print_list[key].append(str(np.nanmax(getattr(self,key).values)))
                except TypeError:
                    print_list[key].append("--")
            maxlen4 = max(maxlen4,len(print_list[key][3]))
            maxlen5 = max(maxlen5,len(print_list[key][4]))
        print("The variables are:")
        print("Name".ljust(maxlen1)+" Type".ljust(maxlen2)+"  Unit".ljust(maxlen3)+"     Min".ljust(maxlen4)+"      Max".ljust(maxlen5))
        for key in sorted(key_list):
            print(print_list[key][0].ljust(maxlen1)+" "+print_list[key][1].ljust(maxlen2)+" ["+print_list[key][2].ljust(maxlen3)+"] "+\
                  print_list[key][3].ljust(maxlen4)+" "+print_list[key][4].ljust(maxlen5))
            #print(key.ljust(maxlen1)+" "+getattr(self,key).kind.ljust(maxlen2)+\
                  #" ["+getattr(self,key).unit.ljust(maxlen3)+"] "+\
                  #str(np.nanmin(getattr(self,key).values)).ljust(maxlen4)+" "+\
                  #str(np.nanmax(getattr(self,key).values)).ljust(maxlen5))
        return
    
    #=======================================================================================
    # The find_center function finds the center in the mesh before loading the full data.
    #=======================================================================================
    def find_center(self,dx,dy,dz):
        
        lc = False
        try: # check if center is defined at all, if not set to (0.5,0.5,0.5)
            lc = len(self.info["center"])
        except TypeError: # No center defined: set to (0.5,0.5,0.5)
            xc = yc = zc = 0.5
        if lc:
            try: # check if center contains numbers
                self.info["center"][0] += 0
                if lc == 3:
                    xc = self.info["center"][0]
                    yc = self.info["center"][1]
                    zc = self.info["center"][2]
                else:
                    print("Bad center format: must have 3 numbers as input.")
                    return
            except TypeError: # if not it should have the format 'sink1', or 'max:density'
                if self.info["center"].startswith("sink"):
                    isink = self.sinks["id"].index(self.info["center"])
                    xc = self.sinks["x"][isink]/self.info["boxlen"]/self.info["unit_l"]
                    yc = self.sinks["y"][isink]/self.info["boxlen"]/self.info["unit_l"]
                    zc = self.sinks["z"][isink]/self.info["boxlen"]/self.info["unit_l"]
                else:
                    xc = yc = zc = 0.5
                    
        return xc,yc,zc
    
    #=======================================================================================
    # The re_center function shifts the coordinates axes around a center. If center="auto"
    # then the function find the cell with the highest density.
    #=======================================================================================
    def re_center(self,newcenter=None):
        
        #try: # check if newcenter is defined
            #lc = len(newcenter)
            #self.data["x"]["values"] = (self.data["x"]["values"] + self.info["xc"])*conf.constants[self.info["scale"]]
            #self.data["y"]["values"] = (self.data["y"]["values"] + self.info["yc"])*conf.constants[self.info["scale"]]
            #if self.info["ndim"] > 2:
                #self.data["z"]["values"] = (self.data["z"]["values"] + self.info["zc"])*conf.constants[self.info["scale"]]
            
            ## Re-scale the cell and box sizes
            #self.data["dx"]["values"] = self.data["dx"]["values"]*conf.constants[self.info["scale"]]
            #self.info["boxsize"] = self.info["boxsize"]*conf.constants[self.info["scale"]]
            
            ## Re-center sinks
            #if self.info["nsinks"] > 0:
                #self.sinks["x"     ] = (self.sinks["x"]+self.info["xc"])*conf.constants[self.info["scale"]]
                #self.sinks["y"     ] = (self.sinks["y"]+self.info["yc"])*conf.constants[self.info["scale"]]
                #self.sinks["z"     ] = (self.sinks["z"]+self.info["zc"])*conf.constants[self.info["scale"]]
                #self.sinks["radius"] =  self.sinks["radius"]/self.info["boxsize"]
            
                ##for key in self.sinks.keys():
                    ##self.sinks[key]["x"     ] = (self.sinks[key]["x"]+self.info["xc"])*conf.constants[self.info["scale"]]
                    ##self.sinks[key]["y"     ] = (self.sinks[key]["y"]+self.info["yc"])*conf.constants[self.info["scale"]]
                    ##self.sinks[key]["z"     ] = (self.sinks[key]["z"]+self.info["zc"])*conf.constants[self.info["scale"]]
                    ##self.sinks[key]["radius"] = self.sinks[key]["radius"]/self.info["boxsize"]
            
            #self.info["center"] = newcenter
        
        try: # check if newcenter is defined
            lc = len(newcenter)
            self.info["center"] = newcenter
        except TypeError:
            pass
        

        try: # check if center is defined at all, if not set to (0.5,0.5,0.5)
            lc = len(self.info["center"])
            try: # check if center contains numbers
                self.info["center"][0] += 0
                if lc == 3:
                    xc = self.info["center"][0]*self.info["boxsize"]
                    yc = self.info["center"][1]*self.info["boxsize"]
                    zc = self.info["center"][2]*self.info["boxsize"]
                else:
                    print("Bad center format: must have 3 numbers as input.")
                    return
            except TypeError: # if not it should have the format 'sink1', or 'max:density'
                if self.info["center"].startswith("sink"):
                    isink = self.sinks["id"].index(self.info["center"])
                    xc = self.sinks["x"][isink]
                    yc = self.sinks["y"][isink]
                    zc = self.sinks["z"][isink]
                elif self.info["center"].startswith("max"):
                    cvar=self.info["center"].split(":")[1]
                    maxloc = np.argmax(getattr(self,cvar).values)
                    xc = self.x_raw.values[maxloc]
                    yc = self.y_raw.values[maxloc]
                    zc = self.z_raw.values[maxloc]
                elif self.info["center"].startswith("min"):
                    cvar=self.info["center"].split(":")[1]
                    minloc = np.argmin(getattr(self,cvar).values)
                    xc = self.x_raw.values[minloc]
                    yc = self.y_raw.values[minloc]
                    zc = self.z_raw.values[minloc]
                elif self.info["center"].startswith("av"):
                    cvar=self.info["center"].split(":")[1]
                    [op_parsed,depth,status] = self.parse_operation(cvar)
                    select = eval("np.where("+op_parsed+")")
                    xc = np.average(self.x_raw.values[select])
                    yc = np.average(self.y_raw.values[select])
                    zc = np.average(self.z_raw.values[select])
                else:
                    print("Bad center value:"+str(self.info["center"]))
                    return
                
        except TypeError: # No center defined: set to (0.5,0.5,0.5)
            xc = yc = zc = 0.5*self.info["boxsize"]

        self.x.values = (self.x_raw.values - xc)/conf.constants[self.info["scale"]]
        if self.info["ndim"] > 1:
            self.y.values = (self.y_raw.values - yc)/conf.constants[self.info["scale"]]
        if self.info["ndim"] > 2:
            self.z.values = (self.z_raw.values - zc)/conf.constants[self.info["scale"]]
        self.info["xc"] = xc/conf.constants[self.info["scale"]]
        self.info["yc"] = yc/conf.constants[self.info["scale"]]
        self.info["zc"] = zc/conf.constants[self.info["scale"]]
        
        # Re-scale the cell and box sizes
        self.dx.values = self.dx_raw.values/conf.constants[self.info["scale"]]
        self.info["boxsize_scaled"] = self.info["boxsize"]/conf.constants[self.info["scale"]]
        
        # Re-center sinks
        if self.info["nsinks"] > 0:
            #for key in self.sinks.keys():
            self.sinks["x"     ] = self.sinks["x"]/conf.constants[self.info["scale"]]-self.info["xc"]
            self.sinks["y"     ] = self.sinks["y"]/conf.constants[self.info["scale"]]-self.info["yc"]
            self.sinks["z"     ] = self.sinks["z"]/conf.constants[self.info["scale"]]-self.info["zc"]
            self.sinks["radius"] = self.sinks["radius"]*self.info["boxsize"]/conf.constants[self.info["scale"]]
        
        return
        
    #=======================================================================================
    # This function reads the sink particle data if present.
    #=======================================================================================
    def read_sinks(self):
        
        sinkfile = self.info["infile"]+"/sink_"+self.info["infile"].split("_")[-1]+".csv"
        try:
            sinklist = np.loadtxt(sinkfile,delimiter=",")
        except IOError:
            self.info["nsinks"] = 0
            return
        if np.shape(sinklist)[0] == 0:
            self.info["nsinks"] = 0
        else:
            list_shape = np.shape(np.shape(sinklist))[0]
            if list_shape == 1:
                sinklist = np.reshape(sinklist, (1, np.shape(sinklist)[0]))
                self.info["nsinks"] = 1
            else:
                self.info["nsinks"] = np.shape(sinklist)[0]
            try:
                r_sink = self.info["ir_cloud"]/(2.0**self.info["levelmax"])
            except KeyError:
                try:
                    r_sink = self.info["ncell_racc"]/(2.0**self.info["levelmax"])
                except KeyError:
                    r_sink = 4.0/(2.0**self.info["levelmax"])
            self.sinks = dict()
            j = 0
            for entry in conf.default_values["sink_format"]:
                self.sinks[entry] = sinklist[:,j]
                j += 1
            self.sinks["x"] *= self.info["unit_l"]
            self.sinks["y"] *= self.info["unit_l"]
            self.sinks["z"] *= self.info["unit_l"]
            self.sinks["radius"] = np.full(self.info["nsinks"],r_sink)
            ids = []
            for i in range(self.info["nsinks"]):
                ids.append("sink"+str(int(self.sinks["number"][i])))
            self.sinks["id"] = ids
            #print("Read %i sink particles" % self.info["nsinks"])
            
        return
            
    #=======================================================================================
    # The update_values function reads in a new ramses output and updates the fields in an
    # existing data structure. It also updates all the derived variables at the same time.
    #=======================================================================================
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
    # The function get_units returns the appropriate scaling for a variable which was read
    # in code units by the data loader. It tries to identify if we are dealing with a
    # density or a pressure and returns the appropriate combination of ud, ul and ut. It
    # also returns the unit as a string for plotting on the axes.
    #=======================================================================================
    def get_units(self,string,ud,ul,ut,scale="cm"):
        if string == "density":
            return [ud,"g/cm3"]
        elif string.startswith("velocity"):
            return [ul/ut,"cm/s"]
        elif string.startswith("momentum"):
            return [ud*ul/ut,"g/cm2/s"]
        elif string.startswith("B_"):
            return [np.sqrt(4.0*np.pi*ud*(ul/ut)**2),"G"]
        elif string == "thermal_pressure":
            return [ud*((ul/ut)**2),"erg/cm3"]
        elif string == "total_energy":
            return [ud*((ul/ut)**2),"erg/cm3"]
        elif string.startswith("radiative_energy"):
            return [ud*((ul/ut)**2),"erg/cm3"]
        elif string == "x":
            return [ul,scale]
        elif string == "y":
            return [ul,scale]
        elif string == "z":
            return [ul,scale]
        elif string == "dx":
            return [ul,scale]
        elif string == "temperature":
            return [1.0,"K"]
        else:
            for key in conf.default_units.keys():
                if string == key:
                    new_string = conf.default_units[string][0].replace("unit_d","self.info[\"unit_d\"]")
                    new_string = new_string.replace("unit_l","self.info[\"unit_l\"]")
                    new_string = new_string.replace("unit_t","self.info[\"unit_t\"]")
                    uu = eval(new_string)
                    return [uu,conf.default_units[string][1]]
            return [1.0,""]


    #=======================================================================================
    # The new field function is used to create a new data field. Say you want to take the
    # log of the density. You create a new field by calling:
    # mydata.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    # The operation string is then evaluated using the 'eval' function.
    #=======================================================================================
    def new_field(self,name,operation="",unit="",label="",verbose=True,values=[],norm=1.0,kind="scalar",\
                  vec_x=False,vec_y=False,vec_z=False,update=False):
        
        # Case where values are given and no operation is to be computed
        if (len(operation) == 0) and (len(values) > 0):
            new_data = values
            op_parsed = operation
            depth = -1
            if hasattr(self,name):
                if verbose:
                    print("Warning: field "+name+" already exists and will be overwritten.")
                theField = getattr(self,name)
                theField.values = values
                if not update:
                    theField.unit = unit
                    theField.label = label
                    theField.operation = operation
                    theField.depth = depth
                    theField.norm = norm
                    theField.kind = kind
                    theField.parent = parent
                    if vec_x:
                        theField.x = vec_x
                    if vec_y:
                        theField.y = vec_y
                    if vec_z:
                        theField.z = vec_z
            else:
                dataField = OsirisData(values=new_data,unit=unit,label=label,operation=op_parsed,depth=depth+1,\
                                       norm=norm,kind=kind,parent=self,vec_x=vec_x,vec_y=vec_y,vec_z=vec_z,name=name)
                setattr(self, name, dataField)
            
        # Case where operation is required
        elif (len(operation) > 0) and (len(values) == 0):
            [op_parsed,depth,status] = self.parse_operation(operation)
            if status == 0:
                print("Cannot combine scalar and vector fields.")
                return
            elif status == 1:
                try:
                    new_data = eval(op_parsed)
                except NameError:
                    if verbose:
                        print("Error parsing operation when trying to create variable: "+name)
                        print("The attempted operation was: "+op_parsed)
                    return
                dataField = OsirisData(values=new_data,unit=unit,label=label,operation=op_parsed,depth=depth+1,\
                               norm=norm,kind=kind,parent=self,name=name)
                if hasattr(self,name) and verbose:
                    print("Warning: field "+name+" already exists and will be overwritten.")
                setattr(self, name, dataField)
            elif status == 2:
                # Dealing with vector fields: first create x,y,z components
                comps = ["_x","_y","_z"]
                for n in range(self.info["ndim"]):
                    [op_parsed,depth,stat_n] = self.parse_operation(operation,suffix=comps[n])
                    if stat_n == 1:
                        try:
                            new_data = eval(op_parsed)
                        except NameError:
                            if verbose:
                                print("Error parsing operation when trying to create variable: "+name+comps[n])
                                print("The attempted operation was: "+op_parsed)
                            return
                        dataField = OsirisData(values=new_data,unit=unit,label=label,operation=op_parsed,depth=depth+1,\
                                       norm=norm,kind=kind,parent=self,name=name)
                        if hasattr(self,name+comps[n]) and verbose:
                            print("Warning: field "+name+comps[n]+" already exists and will be overwritten.")
                        setattr(self, name+comps[n], dataField)
                    else:
                        print("Error: failed to create vector field.")
                        return
                # Dealing with vector fields: then create vector container
                self.vector_field(name=name,key=name)
        
        # Case where both values and operation are empty
        elif (len(operation) == 0) and (len(values) == 0):
            dataField = OsirisData(unit=unit,label=label,parent=self,name=name)
            setattr(self, name, dataField)
        # Case where both values and operation are required
        else:
            print("Both values and operation are defined. Please choose only one.")
        
        return
    
    #=======================================================================================
    # Delete a variable field from the memory
    #=======================================================================================
    def delete_field(self,name):
        
        delattr(self,name)
        
        return
    
    #=======================================================================================
    # The operation parser converts an operation string into an expression which contains
    # variables from the data dictionary. If a name from the variable list, e.g. "density",
    # is found in the operation, it is replaced by self.density.values so that it
    # can be properly evaluated by the 'eval' function in the 'new_field' function.
    #=======================================================================================
    def parse_operation(self,operation,suffix=""):
        
        max_depth = 0
        # Add space before and after to make it easier when searching for characters before
        # and after
        expression = " "+operation+" "
        # Sort the list of variable keys in the order of the longest to the shortest.
        # This guards against replacing 'B' inside 'logB' for example.
        key_list = self.get_var_list()
        key_list = sorted(key_list,key=lambda x:len(x),reverse=True)
        # For replacing, we need to create a list of hash keys to replace on instance at a
        # time
        hashkeys  = dict()
        hashcount = 0
        found_scalar = False
        found_vector = False
        for key in key_list:
            # Search for all instances in string
            loop = True
            loc = 0
            while loop:
                loc = expression.find(key,loc)
                if loc == -1:
                    loop = False
                else:
                    # Check character before and after. If they are either a letter or a '_'
                    # then the instance is actually part of another variable or function
                    # name.
                    char_before = expression[loc-1]
                    char_after  = expression[loc+len(key)]
                    bad_before = (char_before.isalpha() or (char_before == "_"))
                    bad_after = (char_after.isalpha() or (char_after == "_"))
                    hashcount += 1
                    if (not bad_before) and (not bad_after):
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        # Store the data key in the hash table
                        #hashkeys[theHash] = "self.data[\""+key+"\"][\"values\"]"
                        hashkeys[theHash] = "self.get(\""+key+suffix+"\")"
                        expression = expression.replace(key,theHash,1)
                        max_depth = max(max_depth,getattr(self,key+suffix).depth)
                        if getattr(self,key+suffix).kind == "scalar":
                            found_scalar = True
                        if getattr(self,key+suffix).kind == "vector":
                            found_vector = True
                    else:
                        # Replace anyway to prevent from replacing "x" in "max("
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        hashkeys[theHash] = key
                        expression = expression.replace(key,theHash,1)
                    loc += 1
        # Now go through all the hashes in the table and build the final expression
        for theHash in hashkeys.keys():
            expression = expression.replace(theHash,hashkeys[theHash])
        
        if found_scalar and found_vector:
            status = 0
        elif found_scalar:
            status = 1
        elif found_vector:
            status = 2
        else:
            status = 3
        
        return [expression,max_depth,status]
    
    #=======================================================================================
    # The function get returns the values of the selected variable
    #=======================================================================================
    def get(self,var):
        return getattr(getattr(self,var),'values')
    
    #=======================================================================================
    # The function returns the list of variables
    #=======================================================================================
    def get_var_list(self,types=False):
        key_list = []
        typ_list = []
        att_list =  dir(self)
        for att in att_list:
            class_name = getattr(self,att).__class__.__name__
            if class_name == 'OsirisData':
                key_list.append(att)
                typ_list.append(getattr(self,att).kind)
        if types:
            return [key_list,typ_list]
        else:
            return key_list
    
    #=======================================================================================
    # Create a hash table for all the cells in the domain
    #=======================================================================================
    def create_hash_table(self):
        
        print("Building hash table")
        self.hash_table = dict()
        for icell in range(self.info["ncells"]):
            igrid = int(self.x_box.values[icell]/self.dx_box.values[icell])
            jgrid = int(self.y_box.values[icell]/self.dx_box.values[icell])
            kgrid = int(self.z_box.values[icell]/self.dx_box.values[icell])
            theHash = str(igrid)+','+str(jgrid)+','+str(kgrid)+','+str(int(self.level.values[icell]))
            self.hash_table[theHash] = icell

        return

    #=======================================================================================
    # Create dummy variables containing the components of the vectors
    #=======================================================================================
    def create_vector_containers(self):
    
        list_vars = self.get_var_list()
    
        if self.info["ndim"] > 1:
            for i in range(len(list_vars)):
                key = list_vars[i]
                if key.endswith("_x"):
                    rawkey = key[:-2]
                    ok = True
                    try:
                        k1 = len(self.get(rawkey+"_y"))
                    except AttributeError:
                        ok = False
                    if self.info["ndim"] > 2:
                        try:
                            k2 = len(self.get(rawkey+"_z"))
                        except AttributeError:
                            ok = False
                    
                    if ok:
                        vec_name = rawkey
                        while hasattr(self,vec_name):
                            vec_name += "_vec"
                        self.vector_field(name=vec_name,key=rawkey)

        return

    #=======================================================================================
    # Create vector field
    #=======================================================================================
    def vector_field(self,name="",key=""):
    
        v_x=getattr(self,key+"_x")
        v_y=getattr(self,key+"_y")
        v_x.vector_component = True
        v_y.vector_component = True
        v_z = False
        if self.info["ndim"] > 2:
            v_z=getattr(self,key+"_z")
            v_z.vector_component = True
        self.new_field(name=name,values="--",label=name,vec_x=v_x,vec_y=v_y,vec_z=v_z,kind="vector",unit=v_x.unit)
        
        return

#=======================================================================================
#=======================================================================================
# End of class LoadRamsesData()
#=======================================================================================
#=======================================================================================


#=======================================================================================
#=======================================================================================
# USEFUL TOOLS
#=======================================================================================
#=======================================================================================

#=======================================================================================
# Determine binary offset when reading fortran binary files and return unpacked data
#=======================================================================================
def get_binary_data(fmt="",ninteg=0,nlines=0,nfloat=0,nstrin=0,nquadr=0,content=None,correction=0):
    
    offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4 + correction
    byte_size = {"i":4,"d":8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = eval(fmt[0:len(fmt)-1])
    pack_size = mult*byte_size[fmt[-1]]
    
    return struct.unpack(fmt, content[offset:offset+pack_size])
    
