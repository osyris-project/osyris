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
#import plot_osiris
import config_osiris as conf

divider = "============================================"



class OsirisData():
    
    def __init__(self,the_values,the_unit,the_label,the_operation,the_depth):
        
        setattr(self, 'values', the_values)
        setattr(self, 'unit'  , the_unit)
        setattr(self, 'label', the_label)
        setattr(self, 'operation', the_operation)
        setattr(self, 'depth', the_depth)
        
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
        
        ## Read in custom variables if any from the configuration file
        #conf.additional_variables(self)
        
        ## Print exit message
        #print("Memory used: %.2f Mb" % (len(self.data)*self.info["ncells"]*8.0/1.0e6))
        #print(self.info["infile"]+" successfully loaded")
        #if verbose:
            #self.print_info()
        #print(divider)
    
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
        try:
            with open(infofile) as f:
                content = f.readlines()
            f.close()
        except IOError:
            # Clean exit if the file was not found
            print("Info file not found: "+infofile)
            return 0
        
        if not update:
            self.info = dict()
        for line in content:
            sp = line.split("=")
            if len(sp) > 1:
                try:
                    self.info[sp[0].strip()] = eval(sp[1].strip())
                except NameError:
                    self.info[sp[0].strip()] = sp[1].strip()
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
        self.info["nout"     ] = nout
        
        print(divider)
        
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
                content.append("variable #"+str(ivar)+" : "+var)
        # Read the total number of hydro variables
        for line in content:
            sp = line.split("=")
            if len(sp) > 1:
                if sp[0].strip() == "nvar":
                    self.info["nvar"] = int(sp[1].strip())
                    break
        # Now go through all the variables and check if they are to be read or skipped
        var_read = np.ones([self.info["nvar"]+5],dtype=np.bool)
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
        #var_read += "1 1 1 1 1 "
        list_vars.extend(("level","x","y","z","dx"))
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
                                
                # Loop over domains
                for j in range(nboundary+self.info["ncpu"]):
                    
                    ncache = ngridlevel[j,ilevel]
                    
                    # Skip two lines of integers
                    nlines_hydro += 2
                    ninteg_hydro += 2
                    
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
                                for ivar in range(self.info["nvar"]):
                                    if var_read[ivar]:
                                        offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*self.info["nvar"]+ivar)*(ncache+1)) + nstrin_hydro + 4
                                        var[:ncache,ind,jvar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                        jvar += 1
                                var[:ncache,ind,-5] = float(ilevel+1)
                                for n in range(self.info["ndim"]):
                                    xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                    var[:ncache,ind,-4+n] = xyz[:ncache,ind,n]*self.info["boxlen"]
                                var[:ncache,ind,-1] = dxcell*self.info["boxlen"]
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
                        
                        nfloat_hydro += ncache*twotondim*self.info["nvar"]
                        nlines_hydro += twotondim*self.info["nvar"]
                
                # Now increment the offsets while looping through the levels
                ninteg1 = ninteg_amr
                nfloat1 = nfloat_amr
                nlines1 = nlines_amr
                nstrin1 = nstrin_amr
                
                ninteg2 = ninteg_hydro
                nfloat2 = nfloat_hydro
                nlines2 = nlines_hydro
                nstrin2 = nstrin_hydro
        
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
        if not update:
            self.data = dict()
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            if not update:
                self.data[theKey] = dict()
            [norm,uu] = self.get_units(theKey,self.info["unit_d"],self.info["unit_l"],self.info["unit_t"],self.info["scale"])
            # Replace "_" with " " to avoid error with latex when saving figures
            theLabel = theKey.replace("_"," ")
            # Use the 'new_field' function to create data field
            self.new_field(name=theKey,operation="",unit=uu,label=theLabel,values=master_data_array[:,i]*norm,verbose=False)
        
        ## Re-center the mesh around chosen center
        #self.re_center()

        return 1
        
    #=======================================================================================
    # Print information about the data that was loaded.
    #=======================================================================================
    def print_info(self):
        print("--------------------------------------------")
        for key in sorted(self.info.keys()):
            print(key+": "+str(self.info[key]))
        print("--------------------------------------------")
        maxlen1 = 0
        maxlen2 = 0
        maxlen3 = 0
        maxlen4 = 0
        for key in sorted(self.data.keys()):
            maxlen1 = max(maxlen1,len(key))
            maxlen2 = max(maxlen2,len(self.data[key]["unit"]))
            maxlen3 = max(maxlen3,len(str(np.nanmin(self.data[key]["values"]))))
            maxlen4 = max(maxlen4,len(str(np.nanmax(self.data[key]["values"]))))
        print("The variables are:")
        print("Name".ljust(maxlen1)+" "+"Unit".ljust(maxlen2)+"   Min".ljust(maxlen3)+"    Max".ljust(maxlen4))
        for key in sorted(self.data.keys()):
            print(key.ljust(maxlen1)+" ["+self.data[key]["unit"].ljust(maxlen2)+"] "+\
                  str(np.nanmin(self.data[key]["values"])).ljust(maxlen3)+" "+\
                  str(np.nanmax(self.data[key]["values"])).ljust(maxlen4))
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
        
        try: # check if newcenter is defined
            lc = len(newcenter)
            self.data["x"]["values"] = (self.data["x"]["values"] + self.info["xc"])*conf.constants[self.info["scale"]]
            self.data["y"]["values"] = (self.data["y"]["values"] + self.info["yc"])*conf.constants[self.info["scale"]]
            if self.info["ndim"] > 2:
                self.data["z"]["values"] = (self.data["z"]["values"] + self.info["zc"])*conf.constants[self.info["scale"]]
            
            # Re-scale the cell and box sizes
            self.data["dx"]["values"] = self.data["dx"]["values"]*conf.constants[self.info["scale"]]
            self.info["boxsize"] = self.info["boxsize"]*conf.constants[self.info["scale"]]
            
            # Re-center sinks
            if self.info["nsinks"] > 0:
                self.sinks["x"     ] = (self.sinks["x"]+self.info["xc"])*conf.constants[self.info["scale"]]
                self.sinks["y"     ] = (self.sinks["y"]+self.info["yc"])*conf.constants[self.info["scale"]]
                self.sinks["z"     ] = (self.sinks["z"]+self.info["zc"])*conf.constants[self.info["scale"]]
                self.sinks["radius"] =  self.sinks["radius"]/self.info["boxsize"]
            
                #for key in self.sinks.keys():
                    #self.sinks[key]["x"     ] = (self.sinks[key]["x"]+self.info["xc"])*conf.constants[self.info["scale"]]
                    #self.sinks[key]["y"     ] = (self.sinks[key]["y"]+self.info["yc"])*conf.constants[self.info["scale"]]
                    #self.sinks[key]["z"     ] = (self.sinks[key]["z"]+self.info["zc"])*conf.constants[self.info["scale"]]
                    #self.sinks[key]["radius"] = self.sinks[key]["radius"]/self.info["boxsize"]
            
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
                    maxloc = np.argmax(self.data[cvar]["values"])
                    xc = self.data["x"]["values"][maxloc]
                    yc = self.data["y"]["values"][maxloc]
                    zc = self.data["z"]["values"][maxloc]
                elif self.info["center"].startswith("min"):
                    cvar=self.info["center"].split(":")[1]
                    minloc = np.argmin(self.data[cvar]["values"])
                    xc = self.data["x"]["values"][minloc]
                    yc = self.data["y"]["values"][minloc]
                    zc = self.data["z"]["values"][minloc]
                elif self.info["center"].startswith("av"):
                    cvar=self.info["center"].split(":")[1]
                    [op_parsed,depth] = self.parse_operation(cvar)
                    select = eval("np.where("+op_parsed+")")
                    xc = np.average(self.data["x"]["values"][select])
                    yc = np.average(self.data["y"]["values"][select])
                    zc = np.average(self.data["z"]["values"][select])
                else:
                    print("Bad center value:"+str(self.info["center"]))
                    return
                
        except TypeError: # No center defined: set to (0.5,0.5,0.5)
            xc = yc = zc = 0.5*self.info["boxsize"]

        self.data["x"]["values"] = (self.data["x"]["values"] - xc)/conf.constants[self.info["scale"]]
        if self.info["ndim"] > 1:
            self.data["y"]["values"] = (self.data["y"]["values"] - yc)/conf.constants[self.info["scale"]]
        if self.info["ndim"] > 2:
            self.data["z"]["values"] = (self.data["z"]["values"] - zc)/conf.constants[self.info["scale"]]
        self.info["xc"] = xc/conf.constants[self.info["scale"]]
        self.info["yc"] = yc/conf.constants[self.info["scale"]]
        self.info["zc"] = zc/conf.constants[self.info["scale"]]
        
        # Re-scale the cell and box sizes
        self.data["dx"]["values"] = self.data["dx"]["values"]/conf.constants[self.info["scale"]]
        self.info["boxsize"] = self.info["boxsize"]/conf.constants[self.info["scale"]]
        
        # Re-center sinks
        if self.info["nsinks"] > 0:
            #for key in self.sinks.keys():
            self.sinks["x"     ] = self.sinks["x"]/conf.constants[self.info["scale"]]-self.info["xc"]
            self.sinks["y"     ] = self.sinks["y"]/conf.constants[self.info["scale"]]-self.info["yc"]
            self.sinks["z"     ] = self.sinks["z"]/conf.constants[self.info["scale"]]-self.info["zc"]
            self.sinks["radius"] = self.sinks["radius"]*self.info["boxsize"]
        
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
    def update_values(self,nout="none",lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="",\
                      path="",variables=[],verbose=False):
        
        # Check if new output number is requested. If not, use same nout as before
        if nout == "none":
            nout = self.info["nout"]
        
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
        key_list = sorted(self.data.keys(),key=lambda x:self.data[x]["depth"])
        for key in key_list:
            if len(self.data[key]["operation"]) > 0:
                print("Re-computing "+key)
                self.data[key]["values"] = eval(self.data[key]["operation"])
        
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
            return [1.0,""]











    #=======================================================================================
    # The new field function is used to create a new data field. Say you want to take the
    # log of the density. You create a new field by calling:
    # mydata.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    # The operation string is then evaluated using the 'eval' function.
    #=======================================================================================
    def new_field(self,name,operation="",unit="",label="",verbose=True,values=[]):
        
        if (len(operation) == 0) and (len(values) > 0):
            new_data = values
            op_parsed = operation
            depth = -1
        else:
            [op_parsed,depth] = self.parse_operation(operation)
            try:
                new_data = eval(op_parsed)
            except NameError:
                if verbose:
                    print("Error parsing operation when trying to create variable: "+name)
                    print("The attempted operation was: "+op_parsed)
                return
        #TheDict = dict()
        #TheDict["values"   ] = new_data
        #TheDict["unit"     ] = unit
        #TheDict["label"    ] = label
        #TheDict["operation"] = op_parsed
        #TheDict["depth"    ] = depth+1
        
        dataField = OsirisData(new_data,unit,label,op_parsed,depth+1)
        
        
        setattr(self, name, dataField)
        
        return
    
    #=======================================================================================
    # The operation parser converts an operation string into an expression which contains
    # variables from the data dictionary. If a name from the variable list, e.g. "density",
    # is found in the operation, it is replaced by self.data["density"]["values"] so that it
    # can be properly evaluated by the 'eval' function in the 'new_field' function.
    #=======================================================================================
    def parse_operation(self,operation):
        
        max_depth = 0
        # Add space before and after to make it easier when searching for characters before
        # and after
        expression = " "+operation+" "
        # Sort the list of variable keys in the order of the longest to the shortest.
        # This guards against replacing 'B' inside 'logB' for example.
        key_list = sorted(self.data.keys(),key=lambda x:len(x),reverse=True)
        # For replacing, we need to create a list of hash keys to replace on instance at a
        # time
        hashkeys  = dict()
        hashcount = 0
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
                        hashkeys[theHash] = "self.get(\""+key+"\")"
                        expression = expression.replace(key,theHash,1)
                        max_depth = max(max_depth,self.data[key]["depth"])
                    else:
                        # Replace anyway to prevent from replacing "x" in "max("
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        hashkeys[theHash] = key
                        expression = expression.replace(key,theHash,1)
                    loc += 1
        # Now go through all the hashes in the table and build the final expression
        for theHash in hashkeys.keys():
            expression = expression.replace(theHash,hashkeys[theHash])
    
        return [expression,max_depth]
    
    #=======================================================================================
    # The function get returns the values of the selected variable
    #=======================================================================================
    def get(self,var):
        return getattr(getattr(self,var),'values')
    
        
    #=======================================================================================
    # This function writes the RAMSES data to a VTK file for 3D visualization
    #=======================================================================================
    def to_vtk(self,fname="osiris_data.vtu",variables=False):
        
        try:
            from scipy.spatial import Delaunay
        except ImportError:
            print("Scipy Delaunay library not found. This is needed for VTK output. Exiting.")

        # Print status
        if not fname.endswith(".vtu"):
            fname += ".vtu"
        print("Writing data to VTK file: "+fname)
        
        # Coordinates ot RAMSES cell centers
        points = np.array([self.get("x"),self.get("y"),self.get("z")]).T
        
        # Compute Delaunay tetrahedralization from cell nodes
        # Note that this step can take a lot of time!
        ncells = self.info["ncells"]
        print("Computing Delaunay mesh with %i points." % ncells)
        print("This may take some time...")
        tri = Delaunay(points)
        ntetra = np.shape(tri.simplices)[0]
        nverts = ntetra*4
        print("Delaunay mesh with %i tetrahedra complete." % ntetra)

        # Create list of variables by grouping x,y,z components together
        nvarmax = len(self.data.keys())
        n_components = [] #np.zeros([nvarmax],dtype=np.int32)
        varlist = []
        varnames = []
        for key in self.data.keys():
            if key.endswith("_x") or key.endswith("_y") or key.endswith("_z"):
                rawkey = key[:-2]
                try:
                    k = len(self.get(rawkey+"_x"))+len(self.get(rawkey+"_y"))+len(self.get(rawkey+"_z"))
                    ok = True
                    for i in range(np.shape(varlist)[0]):
                        for j in range(n_components[i]):
                            if key == varlist[i][j]:
                                ok = False
                                break
                    if ok:
                        varlist.append([rawkey+"_x",rawkey+"_y",rawkey+"_z"])
                        n_components.append(3)
                        varnames.append(rawkey+"_vec")
                except KeyError:
                    varlist.append([key,"",""])
                    n_components.append(1)
                    varnames.append(key)
            else:
                varlist.append([key,"",""])
                n_components.append(1)
                varnames.append(key)
        
        nvars = len(n_components)

        # Compute byte sizes
        nbytes_xyz   = 3 * ncells * 8
        nbytes_cellc =     nverts * 4
        nbytes_cello =     ntetra * 4
        nbytes_cellt =     ntetra * 4
        nbytes_vars  = np.zeros([nvars],dtype=np.int32)
        for i in range(nvars):
            nbytes_vars[i] = n_components[i] * ncells * 8

        # Compute byte offsets
        offsets = np.zeros([nvars+4],dtype=np.int64)
        offsets[0] = 0                             # xyz coordinates
        offsets[1] = offsets[0] + 4 + nbytes_xyz   # cell connectivity
        offsets[2] = offsets[1] + 4 + nbytes_cellc # cell offsets
        offsets[3] = offsets[2] + 4 + nbytes_cello # cell types
        offsets[4] = offsets[3] + 4 + nbytes_cellt # first hydro variable
        for i in range(nvars-1):
            offsets[i+5] = offsets[i+4] + 4 + nbytes_vars[i]
            
        # Open file for binary output
        f = open(fname, "wb")

        # Write VTK file header
        f.write('<?xml version=\"1.0\"?>\n')
        f.write('<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n')
        f.write('   <UnstructuredGrid>\n')
        f.write('   <Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n' % (ncells,ntetra))
        f.write('      <Points>\n')
        f.write('         <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%i\" />\n' % offsets[0])
        f.write('      </Points>\n')
        f.write('      <Cells>\n')
        f.write('         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%i\" />\n' % offsets[1])
        f.write('         <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%i\" />\n' % offsets[2])
        f.write('         <DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%i\" />\n' % offsets[3])
        f.write('      </Cells>\n')
        f.write('      <PointData>\n')
        for i in range(nvars):
            f.write('         <DataArray type=\"Float64\" Name=\"'+varnames[i]+'\" NumberOfComponents=\"%i\" format=\"appended\" offset=\"%i\" />\n' % (n_components[i],offsets[i+4]))
        f.write('      </PointData>\n')
        f.write('   </Piece>\n')
        f.write('   </UnstructuredGrid>\n')
        f.write('   <AppendedData encoding=\"raw\">\n')
        f.write('_')

        # Now write data in binary. Every data field is preceded by its byte size.
        
        # x,y,z coordinates of the points
        f.write(struct.pack('<i', *[nbytes_xyz]))
        f.write(struct.pack('<%id'%(ncells*3), *np.ravel(points)))

        # Cell connectivity
        f.write(struct.pack('<i', *[nbytes_cellc]))
        f.write(struct.pack('<%ii'%nverts, *np.ravel(tri.simplices)))

        # Cell offsets
        f.write(struct.pack('<i', *[nbytes_cello]))
        f.write(struct.pack('<%ii'%ntetra, *range(4,ntetra*4+1,4)))

        # Cell types: number 10 is tetrahedron in VTK file format
        f.write(struct.pack('<i', *[nbytes_cellt]))
        f.write(struct.pack('<%ii'%ntetra, *np.full(ntetra, 10,dtype=np.int32)))

        # Hydro variables
        #ivar = 0
        for i in range(nvars):
        #for key in self.data.keys():
            if n_components[i] == 3:
                celldata = np.ravel(np.array([self.get(varlist[i][0]),self.get(varlist[i][1]),self.get(varlist[i][2])]).T)
            else:
                celldata = self.get(varlist[i][0])
            f.write(struct.pack('<i', *[nbytes_vars[i]]))
            f.write(struct.pack('<%id'%(ncells*n_components[i]), *celldata))

        # Close file
        f.write('   </AppendedData>\n')
        f.write('</VTKFile>\n')
        f.close()
        
        # File size
        fsize_raw = offsets[nvars+3] + nbytes_vars[nvars-1]
        if fsize_raw > 1000000000:
            fsize = float(fsize_raw)/1.0e9
            funit = "Gb"
        elif fsize_raw > 1000000:
            fsize = float(fsize_raw)/1.0e6
            funit = "Mb"
        elif fsize_raw > 1000:
            fsize = float(fsize_raw)/1.0e3
            funit = "kb"
        else:
            fsize = float(fsize_raw)
            funit = "b"

        print("File "+fname+(" of size %.1f"%fsize)+funit+" succesfully written.")

        return










