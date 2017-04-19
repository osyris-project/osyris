import numpy as np
import glob
import plot_osiris
import read_ramses_data as rd
import config_osiris as conf

divider = "============================================"

#=======================================================================================
# This is the class which will hold the data that you read from the Ramses output
# It calls "rd.ramses_data" which is a fortran file reader.
# It then stores the data in a dictionary named "data"
#=======================================================================================
class LoadRamsesData(plot_osiris.OsirisData):
 
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
    def __init__(self,nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale=False,verbose=False,path="",variables=[]):
                
        # Load the Ramses data using the loader function
        status = self.data_loader(nout=nout,lmax=lmax,center=center,dx=dx,dy=dy,dz=dz,scale=scale,path=path,variables=variables)
        
        if status == 0:
            return
        
        # Re-center the mesh around chosen center
        self.re_center()
        
        # Read in custom variables if any from the configuration file
        conf.additional_variables(self)
        
        # Print exit message
        print(self.info["infile"]+" successfully loaded")
        if verbose:
            self.print_info()
        print(divider)
    
    #=======================================================================================
    # Generate the file name
    #=======================================================================================
    def generate_fname(self,nout,path=""):
        
        if len(path) > 0:
            if path[-1] != "/":
                path=path+"/"
        
        if nout == -1:
            filelist = sorted(glob.glob(path+"output*"))
            infile = filelist[-1]
        else:
            infile = path+"output_"+str(nout).zfill(5)
            
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
        var_read = ""
        list_vars = []
        for line in content:
            sp = line.split(":")
            if len(sp) > 1:
                if (len(variables) == 0) or (variables.count(sp[1].strip()) > 0):
                    var_read += "1 "
                    list_vars.append(sp[1].strip())
                else:
                    var_read += "0 "
        # Make sure we always read the coordinates
        var_read += "1 1 1 1 1 "
        list_vars.extend(("level","x","y","z","dx"))
        nvar_read = len(list_vars)
        
        # Load sink particles if any
        self.read_sinks()
        
        # Find the center
        xc,yc,zc = self.find_center(dx,dy,dz)
                
        print(divider)
        
        # This calls the Fortran data reader and returns the values into the data1
        # array. First a quick scan is made to count how much memory is needed.
        # It then tries to read a hydro_file_descriptor.txt to get a list of
        # variables. If that file is not found, it makes an educated guess for the file
        # contents.
        [active_lmax,nmaxcells,fail] = rd.quick_amr_scan(infile,lmax)
        if fail: # Clean exit if the file was not found
            print(divider)
            return 0
        [data1,nn,fail] = rd.ramses_data(infile,nmaxcells,nvar_read,lmax,var_read,xc,yc,zc,dx,dy,dz,conf.constants[scale],False)
        if fail: # Clean exit if the file was not found
            print(divider)
            return 0
        
        print("Generating data structure... please wait")
        
        # Store the number of cells
        self.info["ncells"] = nn
        
        # This is the master data dictionary. For each entry, the dict has 5 fields.
        # It loops through the list of variables that it got from the file loader.
        if not update:
            self.data = dict()
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            if not update:
                self.data[theKey] = dict()
            [norm,uu] = self.get_units(theKey,self.info["unit_d"],self.info["unit_l"],self.info["unit_t"],self.info["scale"])
            # Use the 'new_field' function to create data field
            self.new_field(name=theKey,operation="",unit=uu,label=theKey,values=data1[:nn,i]*norm)

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
            maxlen3 = max(maxlen3,len(str(np.amin(self.data[key]["values"]))))
            maxlen4 = max(maxlen4,len(str(np.amax(self.data[key]["values"]))))
        print("The variables are:")
        print("Name".ljust(maxlen1)+" "+"Unit".ljust(maxlen2)+"   Min".ljust(maxlen3)+"    Max".ljust(maxlen4))
        for key in sorted(self.data.keys()):
            print(key.ljust(maxlen1)+" ["+self.data[key]["unit"].ljust(maxlen2)+"] "+\
                  str(np.amin(self.data[key]["values"])).ljust(maxlen3)+" "+\
                  str(np.amax(self.data[key]["values"])).ljust(maxlen4))
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
                    xc = self.sinks[self.info["center"]]["x"]/self.info["boxlen"]/self.info["unit_l"]
                    yc = self.sinks[self.info["center"]]["y"]/self.info["boxlen"]/self.info["unit_l"]
                    zc = self.sinks[self.info["center"]]["z"]/self.info["boxlen"]/self.info["unit_l"]
                else:
                    xc = yc = zc = 0.5
                    #if dx+dy+dz > 0.0:
                        #active_lmax,failed = rd.quick_amr_scan(self.info["infile"])
                        #coarse_lmax = int(0.3*(active_lmax - self.info["levelmin"]) + self.info["levelmin"])
                        #[data1,names,nn,fail] = rd.ramses_data(self.info["infile"],coarse_lmax,xc,yc,zc,0.0,0.0,0.0,conf.constants[self.info["scale"]],True)
                        #temp = dict()
                        #list_vars = names.decode().split()
                        #for i in range(len(list_vars)):
                            #theKey = list_vars[i]
                            #temp[theKey] = data1[:nn,i]
                        #if self.info["center"].startswith("max"):
                            #cvar=self.info["center"].split(":")[1]
                            #maxloc = np.argmax(temp[cvar])
                            #xc = temp["x"][maxloc]/self.info["boxlen"]
                            #yc = temp["y"][maxloc]/self.info["boxlen"]
                            #zc = temp["z"][maxloc]/self.info["boxlen"]
                        #elif self.info["center"].startswith("min"):
                            #cvar=self.info["center"].split(":")[1]
                            #minloc = np.argmin(temp[cvar])
                            #xc = temp["x"][minloc]/self.info["boxlen"]
                            #yc = temp["y"][minloc]/self.info["boxlen"]
                            #zc = temp["z"][minloc]/self.info["boxlen"]
                        #else:
                            #print("Bad center value:"+str(self.info["center"]))
                            #return
        return xc,yc,zc
    
    #=======================================================================================
    # The re_center function shifts the coordinates axes around a center. If center="auto"
    # then the function find the cell with the highest density.
    #=======================================================================================
    def re_center(self):

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
                    xc = self.sinks[self.info["center"]]["x"]
                    yc = self.sinks[self.info["center"]]["y"]
                    zc = self.sinks[self.info["center"]]["z"]
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
            for i in range(self.info["nsinks"]):
                key = "sink"+str(i+1)
                self.sinks[key]["x"     ] = self.sinks[key]["x"]/conf.constants[self.info["scale"]]-self.info["xc"]
                self.sinks[key]["y"     ] = self.sinks[key]["y"]/conf.constants[self.info["scale"]]-self.info["yc"]
                self.sinks[key]["z"     ] = self.sinks[key]["z"]/conf.constants[self.info["scale"]]-self.info["zc"]
                self.sinks[key]["radius"] = self.sinks[key]["radius"]*self.info["boxsize"]
        
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
                sinklist = [sinklist[:],sinklist[:]]
            
            self.info["nsinks"] = int(sinklist[-1][0])
            r_sink = self.info["ir_cloud"]/(2.0**self.info["levelmax"])
            
            self.sinks = dict()
            for i in range(self.info["nsinks"]):
                key = "sink"+str(i+1)
                self.sinks[key] = dict()
                self.sinks[key]["mass"    ] = sinklist[i][ 1]
                self.sinks[key]["dmf"     ] = sinklist[i][ 2]
                self.sinks[key]["x"       ] = sinklist[i][ 3]*self.info["unit_l"]
                self.sinks[key]["y"       ] = sinklist[i][ 4]*self.info["unit_l"]
                self.sinks[key]["z"       ] = sinklist[i][ 5]*self.info["unit_l"]
                self.sinks[key]["vx"      ] = sinklist[i][ 6]
                self.sinks[key]["vy"      ] = sinklist[i][ 7]
                self.sinks[key]["vz"      ] = sinklist[i][ 8]
                self.sinks[key]["period"  ] = sinklist[i][ 9]
                self.sinks[key]["lx"      ] = sinklist[i][10]
                self.sinks[key]["ly"      ] = sinklist[i][11]
                self.sinks[key]["lz"      ] = sinklist[i][12]
                self.sinks[key]["acc_rate"] = sinklist[i][13]
                self.sinks[key]["acc_lum" ] = sinklist[i][14]
                self.sinks[key]["age"     ] = sinklist[i][15]
                self.sinks[key]["int_lum" ] = sinklist[i][16]
                self.sinks[key]["Teff"    ] = sinklist[i][17]
                self.sinks[key]["radius"  ] = r_sink
            
            #print("Read %i sink particles" % self.info["nsinks"])
            
        return
            
    #=======================================================================================
    # The new field function is used to create a new data field. Say you want to take the
    # log of the density. You create a new field by calling:
    # mydata.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    # The operation string is then evaluated using the 'eval' function.
    #=======================================================================================
    def new_field(self,name,operation="",unit="",label="",verbose=False,values=[]):
        
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
        self.data[name] = dict()
        self.data[name]["values"   ] = new_data
        self.data[name]["unit"     ] = unit
        self.data[name]["label"    ] = label
        self.data[name]["operation"] = op_parsed
        self.data[name]["depth"    ] = depth+1
        
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
                        hashkeys[theHash] = "self.data[\""+key+"\"][\"values\"]"
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
        
        # Re-center the mesh around chosen center
        self.re_center()
        
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
    # The function get returns the values of the selected variable
    #=======================================================================================
    def get(self,var):
        return self.data[var]["values"]
