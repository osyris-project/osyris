# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np
import struct
from . import config as conf

#=======================================================================================
# This is the class which will holds a scalar or vector field.
#=======================================================================================
class OsyrisField():

    def __init__(self,values=None,unit=None,label=None,operation=None,depth=None,norm=1.0,\
                 parent=None,kind="scalar",vec_x=False,vec_y=False,vec_z=False,name="",\
                 vector_component=False,group="hydro"):

        self.values = values
        self.unit = unit
        self.label = label
        self.operation = operation
        self.depth = depth
        self.norm = norm
        self.kind = kind
        self.parent = parent
        self.name = name
        self.group = group
        if vec_x:
            self.x = vec_x
        if vec_y:
            self.y = vec_y
        if vec_z:
            self.z = vec_z
        self.vector_component = vector_component

        return

    def get(self, only_leafs=True):
        return self.parent.get(self.name, only_leafs=only_leafs)

#=======================================================================================
# This is a dummy class which gives access to common functions to the other
# classes through inheritance.
#=======================================================================================
class OsyrisData:

    def __init__(self):

        return

    #=======================================================================================
    # Read a test file containing information and parameters
    #=======================================================================================
    def read_parameter_file(self,fname="",dict_name="",evaluate=True,verbose=False,delimiter="="):

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
            sp = line.split(delimiter)
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

        # First get maximum length
        maxlen1 = maxlen2 = maxlen3 = maxlen4 = maxlen5 = maxlen6 = 0
        print_list = dict()
        for key in sorted(self.get_var_list()):
            if not getattr(self,key).vector_component:
                print_list[key] = []
                print_list[key].append(key)
                maxlen1 = max(maxlen1,len(key))
                print_list[key].append(getattr(self,key).kind)
                maxlen2 = max(maxlen2,len(print_list[key][1]))
                print_list[key].append(getattr(self,key).group)
                maxlen3 = max(maxlen3,len(print_list[key][2]))
                print_list[key].append(getattr(self,key).unit)
                maxlen4 = max(maxlen4,len(print_list[key][3]))
                print_list[key].append(str(np.nanmin(self.get(key))))
                print_list[key].append(str(np.nanmax(self.get(key))))
                maxlen5 = max(maxlen5,len(print_list[key][4]))
                maxlen6 = max(maxlen6,len(print_list[key][5]))

        # Now print to screen
        rule = "-" * (maxlen1+maxlen2+maxlen3+maxlen4+maxlen5+maxlen6+7)
        print(rule)
        for key in sorted(self.info.keys()):
            theShape = np.shape(self.info[key])
            if len(theShape) > 0:
                try:
                    print(key+": ["+str(self.info[key][0])+" ... "+str(self.info[key][-1])+"]")
                except IndexError:
                    print(key+": "+str(self.info[key]))
            else:
                print(key+": "+str(self.info[key]))
        print(rule)
        print("The variables are:")
        print("Name".ljust(maxlen1)+" Type".ljust(maxlen2)+"  Group".ljust(maxlen3)+\
              " Unit".ljust(maxlen4)+"    Min".ljust(maxlen5)+"     Max".ljust(maxlen6))
        for key in sorted(print_list.keys()):
            print(print_list[key][0].ljust(maxlen1)+" "+print_list[key][1].ljust(maxlen2)+" "+\
                  print_list[key][2].ljust(maxlen3)+" ["+print_list[key][3].ljust(maxlen4)+"] "+\
                  print_list[key][4].ljust(maxlen5)+" "+print_list[key][5].ljust(maxlen6))
        #print(rule)

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
            except TypeError: # if not it should have the format 'sink:1', or 'max:density'
                if self.info["center"].startswith("sink"):
                    isink = np.where(self.sinks["id"] == int(self.info["center"].split(":")[1]))[0][0]
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

            # Find current center
            xc = self.info["xc"] * conf.constants[self.info["scale"]]
            yc = self.info["yc"] * conf.constants[self.info["scale"]]
            zc = self.info["zc"] * conf.constants[self.info["scale"]]

            # Rescale the coordinates
            self.x.values = self.x.values*conf.constants[self.info["scale"]] + xc
            if self.info["ndim"] > 1:
                self.y.values = self.y.values*conf.constants[self.info["scale"]] + yc
            if self.info["ndim"] > 2:
                self.z.values = self.z.values*conf.constants[self.info["scale"]] + zc

            # Re-scale the cell sizes
            self.dx.values = self.dx.values*conf.constants[self.info["scale"]]

            # Re-center sinks
            if self.info["nsinks"] > 0:
                self.sinks["x"     ] = (self.sinks["x"]+self.info["xc"])*conf.constants[self.info["scale"]]
                self.sinks["y"     ] = (self.sinks["y"]+self.info["yc"])*conf.constants[self.info["scale"]]
                self.sinks["z"     ] = (self.sinks["z"]+self.info["zc"])*conf.constants[self.info["scale"]]
                self.sinks["radius"] =  self.sinks["radius"]/self.info["boxsize"]

            # Re-center particles
            if self.info["npart_tot"] > 0:
                self.part_position_x.values = self.part_position_x.values*conf.constants[self.info["scale"]] + xc
                if self.info["ndim"] > 1:
                    self.part_position_y.values = self.part_position_y.values*conf.constants[self.info["scale"]] + yc
                if self.info["ndim"] > 2:
                    self.part_position_z.values = self.part_position_z.values*conf.constants[self.info["scale"]] + zc

            # Store new center in info
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
            except TypeError: # if not it should have the format 'sink:1', or 'max:density'
                cvar=self.info["center"].split(":")[1]
                if self.info["center"].startswith("sink"):
                    isink = np.where(self.sinks["id"] == int(cvar))[0][0]
                    xc = self.sinks["x"][isink]
                    yc = self.sinks["y"][isink]
                    zc = self.sinks["z"][isink]
                elif self.info["center"].startswith("max"):
                    # cvar=self.info["center"].split(":")[1]
                    maxloc = np.argmax(self.get(cvar))
                    xc = self.get("x")[maxloc]
                    yc = self.get("y")[maxloc]
                    zc = self.get("z")[maxloc]
                elif self.info["center"].startswith("min"):
                    # cvar=self.info["center"].split(":")[1]
                    minloc = np.argmin(self.get(cvar))
                    xc = self.get("x")[minloc]
                    yc = self.get("y")[minloc]
                    zc = self.get("z")[minloc]
                elif self.info["center"].startswith("av"):
                    # cvar=self.info["center"].split(":")[1]
                    [op_parsed,depth,grp,status] = self.parse_operation(cvar,only_leafs=True)
                    select = eval("np.where("+op_parsed+")")
                    xc = np.average(self.get("x")[select])
                    yc = np.average(self.get("y")[select])
                    zc = np.average(self.get("z")[select])
                else:
                    print("Bad center value:"+str(self.info["center"]))
                    return

        except TypeError: # No center defined: set to (0.5,0.5,0.5)
            xc = yc = zc = 0.5*self.info["boxsize"]

        self.x.values = (self.x.values - xc)/conf.constants[self.info["scale"]]
        if self.info["ndim"] > 1:
            self.y.values = (self.y.values - yc)/conf.constants[self.info["scale"]]
        if self.info["ndim"] > 2:
            self.z.values = (self.z.values - zc)/conf.constants[self.info["scale"]]
        self.info["xc"] = xc/conf.constants[self.info["scale"]]
        self.info["yc"] = yc/conf.constants[self.info["scale"]]
        self.info["zc"] = zc/conf.constants[self.info["scale"]]

        # Re-scale the cell and box sizes
        self.dx.values = self.dx.values/conf.constants[self.info["scale"]]
        self.info["boxsize_scaled"] = self.info["boxsize"]/conf.constants[self.info["scale"]]

        # Re-center sinks
        if self.info["nsinks"] > 0:
            self.sinks["x"     ] = self.sinks["x"]/conf.constants[self.info["scale"]]-self.info["xc"]
            self.sinks["y"     ] = self.sinks["y"]/conf.constants[self.info["scale"]]-self.info["yc"]
            self.sinks["z"     ] = self.sinks["z"]/conf.constants[self.info["scale"]]-self.info["zc"]
            self.sinks["radius"] = self.sinks["radius"]*self.info["boxsize"]/conf.constants[self.info["scale"]]

        # Re-center particles
        if self.info["npart_tot"] > 0:
            self.part_position_x.values = (self.part_position_x.values - xc)/conf.constants[self.info["scale"]]
            if self.info["ndim"] > 1:
                self.part_position_y.values = (self.part_position_y.values - yc)/conf.constants[self.info["scale"]]
            if self.info["ndim"] > 2:
                self.part_position_z.values = (self.part_position_z.values - zc)/conf.constants[self.info["scale"]]

        return

    #=======================================================================================
    # The new field function is used to create a new data field. Say you want to take the
    # log of the density. You create a new field by calling:
    # mydata.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    # The operation string is then evaluated using the 'eval' function.
    #=======================================================================================
    def new_field(self,name,operation="",unit="",label="",verbose=True,values=[],norm=1.0,kind="scalar",\
                  vec_x=False,vec_y=False,vec_z=False,update=False,group=""):

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
                    theField.parent = self
                    if vec_x:
                        theField.x = vec_x
                    if vec_y:
                        theField.y = vec_y
                    if vec_z:
                        theField.z = vec_z
            else:
                if len(group) == 0:
                    group = "hydro"
                dataField = OsyrisField(values=new_data,unit=unit,label=label,operation=op_parsed,depth=depth+1,\
                                       norm=norm,kind=kind,parent=self,vec_x=vec_x,vec_y=vec_y,vec_z=vec_z,name=name,group=group)
                setattr(self, name, dataField)

        # Case where operation is required
        elif (len(operation) > 0) and (len(values) == 0):
            [op_parsed,depth,grp,status] = self.parse_operation(operation)
            if len(group) == 0:
                    group = grp
            if status == 2: # Only scalar fields
                try:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        new_data = eval(op_parsed)
                except NameError:
                    if verbose:
                        print("Error parsing operation when trying to create variable: "+name)
                        print("The attempted operation was: "+op_parsed)
                    return
                dataField = OsyrisField(values=new_data,unit=unit,label=label,operation=op_parsed,depth=depth+1,\
                               norm=norm,kind=kind,parent=self,name=name,group=group)
                if hasattr(self,name) and verbose:
                    print("Warning: field "+name+" already exists and will be overwritten.")
                setattr(self, name, dataField)
            elif status == 1: # Dealing with vector fields
                # Dealing with vector fields: first create x,y,z components
                comps = ["_x","_y","_z"]
                for n in range(self.info["ndim"]):
                    [op_parsed,depth,grp,stat_n] = self.parse_operation(operation,suffix=comps[n])
                    if stat_n == 2:
                        try:
                            with np.errstate(divide="ignore", invalid="ignore"):
                                new_data = eval(op_parsed)
                        except NameError:
                            if verbose:
                                print("Error parsing operation when trying to create variable: "+name+comps[n])
                                print("The attempted operation was: "+op_parsed)
                            return
                        dataField = OsyrisField(values=new_data,unit=unit,label=label,operation=op_parsed,depth=depth+1,\
                                       norm=norm,kind=kind,parent=self,name=name+comps[n],group=group)
                        if hasattr(self,name+comps[n]) and verbose:
                            print("Warning: field "+name+comps[n]+" already exists and will be overwritten.")
                        setattr(self, name+comps[n], dataField)
                    else:
                        print("Error: failed to create vector field.")
                        return
                # Dealing with vector fields: then create vector container
                self.vector_field(name=name,label=label)

        # Case where both values and operation are empty
        elif (len(operation) == 0) and (len(values) == 0):
            dataField = OsyrisField(unit=unit,label=label,parent=self,name=name,group=group)
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
    # is found in the operation, it is replaced by self.get("density") so that it
    # can be properly evaluated by the 'eval' function in the 'new_field' function.
    #=======================================================================================
    def parse_operation(self,operation,suffix="",only_leafs=False):

        max_depth = 0
        # Add space before and after to make it easier when searching for characters before
        # and after
        expression = " "+operation+" "
        # Sort the list of variable keys in the order of the longest to the shortest.
        # This guards against replacing 'B' inside 'logB' for example.
        key_list = self.get_var_list()
        key_list = sorted(key_list,key=lambda x:len(x),reverse=True)
        # For replacing, we need to create a list of hash keys to replace on instance at a time
        hashkeys  = dict()
        hashcount = 0
        types_found = {"scalar":False,"vector":False,"hydro":False,"amr":False,"grav":False}

        for key in key_list:

            # First look if there are any ".values" in the operation, i.e. vector magnitudes
            keyVal = key+".values"
            if expression.count(keyVal) > 0:
                hashcount += 1
                theHash = "#"+str(hashcount).zfill(5)+"#"
                hashkeys[theHash] = "self."+keyVal
                expression = expression.replace(keyVal,theHash)
                max_depth = max(max_depth,getattr(self,key).depth)
                types_found["scalar"] = True
                types_found[getattr(self,key).group] = True

            # Now search for all instances of individual variables in string
            loop = True
            loc = 0
            while loop:
                loc = expression.find(key,loc)
                if loc == -1:
                    loop = False
                else:
                    # Check character before and after. If they are either a letter or a '_'
                    # then the instance is actually part of another variable or function name.
                    char_before = expression[loc-1]
                    char_after  = expression[loc+len(key)]
                    bad_before = (char_before.isalpha() or (char_before == "_"))
                    bad_after = (char_after.isalpha() or (char_after == "_"))
                    hashcount += 1
                    if (not bad_before) and (not bad_after):
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        # Store the data key in the hash table:
                        if getattr(self,key).kind == "vector":
                            thisKey = key+suffix
                        else:
                            thisKey = key
                        hashkeys[theHash] = "self.get(\""+thisKey+"\",only_leafs="+str(only_leafs)+")"
                        expression = expression.replace(key,theHash,1)
                        max_depth = max(max_depth,getattr(self,thisKey).depth)
                        types_found[getattr(self,thisKey).kind] = True
                        types_found[getattr(self,thisKey).group] = True
                    else:
                        # Replace anyway to prevent from replacing "x" in "max("
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        hashkeys[theHash] = key
                        expression = expression.replace(key,theHash,1)
                    loc += 1
        # Now go through all the hashes in the table and build the final expression
        for theHash in hashkeys.keys():
            expression = expression.replace(theHash,hashkeys[theHash])

        # Determine output group
        if types_found["hydro"]:
            group = "hydro"
        elif types_found["grav"]:
            group = "grav"
        elif types_found["amr"]:
            group = "amr"
        else:
            group = "hydro"

        # Determine exit status
        if types_found["vector"]:
            status = 1
        elif types_found["scalar"]:
            status = 2
        else:
            status = 3

        return [expression,max_depth,group,status]

    #=======================================================================================
    # The function get returns the values of the selected variable.
    # By default, it will only return the leaf cells, but you can choose to return
    # all the cells in the tree by using the argument only_leafs=False.
    #=======================================================================================
    def get(self,var,only_leafs=True):

        # Make sure that we don't use the "only_leafs" indices if we are trying to access
        # particle fields
        if only_leafs and (getattr(self,var).group != "part"):
            return getattr(self,var).values[self.info["leafs"]]
        else:
            return getattr(self,var).values

    #=======================================================================================
    # The function returns the list of variables
    #=======================================================================================
    def get_var_list(self,types=False):
        key_list = []
        typ_list = []
        att_list =  dir(self)
        for att in att_list:
            class_name = getattr(self,att).__class__.__name__
            if class_name == 'OsyrisField':
                key_list.append(att)
                typ_list.append(getattr(self,att).kind)
        if types:
            return [key_list,typ_list]
        else:
            return key_list

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
                        self.vector_field(name=rawkey,label=rawkey)

        return

    #=======================================================================================
    # Create vector field
    #=======================================================================================
    def vector_field(self,name="",values_x=None,values_y=None,values_z=None,unit="",label=""):

        if len(np.shape(values_x)) > 0:
            self.new_field(name+"_x",values=values_x,unit=unit,label=name+"_x",verbose=False)
        if len(np.shape(values_y)) > 0:
            self.new_field(name+"_y",values=values_y,unit=unit,label=name+"_y",verbose=False)
        if len(np.shape(values_z)) > 0:
            self.new_field(name+"_z",values=values_z,unit=unit,label=name+"_z",verbose=False)

        v_x=getattr(self,name+"_x")
        v_y=getattr(self,name+"_y")
        v_x.vector_component = True
        v_y.vector_component = True

        if self.info["ndim"] > 2:
            v_z=getattr(self,name+"_z")
            v_z.vector_component = True
            vals = np.linalg.norm([self.get(name+"_x",only_leafs=False),\
                                   self.get(name+"_y",only_leafs=False),\
                                   self.get(name+"_z",only_leafs=False)],axis=0)
        else:
            v_z = False
            vals = np.linalg.norm([self.get(name+"_x",only_leafs=False),self.get(name+"_y",only_leafs=False)],axis=0)

        self.new_field(name=name,values=vals,label=label,vec_x=v_x,vec_y=v_y,vec_z=v_z,kind="vector",unit=v_x.unit,group=v_x.group)

        return

    #=======================================================================================
    # Get cylindrical basis
    #=======================================================================================
    def get_cylindrical_basis(self,direction):

        pos = np.vstack((self.x.values,self.y.values,self.z.values)).T
        ez   = direction/np.linalg.norm(direction)
        ephi = np.cross(ez,pos)
        ephi_norm = np.linalg.norm(ephi,axis=1)
        ephi = np.vstack((ephi[:,0]/ephi_norm,ephi[:,1]/ephi_norm,ephi[:,2]/ephi_norm)).T
        er   = np.cross(ephi,ez)

        del pos
        return er,ephi,ez

    #=======================================================================================
    # Get cylindrical components of the vector field variable
    #=======================================================================================
    def get_cylindrical_components(self,variable,direction):

        if getattr(self,variable).kind != 'vector':
            print("get_cylindrical_components must be applied to a vector field!")
            return

        if hasattr(self,variable+"_cyl_r"):
            print("****** Warning ****** : cylindrical components of "+variable+" already exist")
            print("you should check if the vector direction was the good one and/or delete the previously computed components")
            return

        er,ephi,ez = self.get_cylindrical_basis(direction)
        vec = np.vstack((getattr(self,variable+"_x").values,getattr(self,variable+"_y").values,getattr(self,variable+"_z").values)).T

        self.new_field(name=variable+"_cyl_r"  ,values=np.sum(vec*er  ,axis=1),unit=getattr(self,variable).unit,label=getattr(self,variable).label+"_r")
        self.new_field(name=variable+"_cyl_phi",values=np.sum(vec*ephi,axis=1),unit=getattr(self,variable).unit,label=getattr(self,variable).label+"_phi")
        self.new_field(name=variable+"_cyl_z"  ,values=np.sum(vec*ez  ,axis=1),unit=getattr(self,variable).unit,label=getattr(self,variable).label+"_z")

        del er,ephi,ez,vec
        return

#=======================================================================================
#=======================================================================================
# End of class OsyrisData()
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
def get_binary_data(fmt="",ninteg=0,nlines=0,nfloat=0,nstrin=0,nquadr=0,nlongi=0,offset=None,content=None,correction=0):

    if offset is None:
        offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16
    offset += 4 + correction
    byte_size = {"b": 1 , "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = eval(fmt[0:len(fmt)-1])
    pack_size = mult*byte_size[fmt[-1]]

    return struct.unpack(fmt, content[offset:offset+pack_size])
