import numpy as np
import glob
import read_ramses_data as rd
import matplotlib.pyplot as plt

#=======================================================================================
# Common variables
#=======================================================================================
scalelist = {"cm": 1.0, "au": 1.495980e+13, "pc": 3.085678e+18}
divider = "============================================"

#=======================================================================================
# This is the class which will hold the data that you read from the Ramses output
# It calls "rd.ramses_data" which is a fortran file reader.
# It then stores the data in a dictionary named "data"
#=======================================================================================
class RamsesOutput:
 
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
    def __init__(self,nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="cm"):
        
        # Generate filename from output number
        if nout == -1:
            filelist = sorted(glob.glob("output*"))
            infile = filelist[-1]
        else:
            infile = "output_"+str(nout).zfill(5)
        
        # Define a center for reading in the data. If it is set to "auto", then it is
        # first set to [0.5,0.5,0.5] for the output reader.
        try:
            center += 0
            xc,yc,zc = center[0:3]
        except TypeError:
            xc = yc = zc = 0.5
        
        print divider
        
        # This calls the Fortran data reader and returns the values into the data1 array.
        # It tries to read a hydro_file_descriptor.txt to get a list of variables.
        # If that file is not found, it makes an educated guess for the file contents.
        [data1,names,nn,ncpu,ndim,levelmin,levelmax,nstep,boxsize,time,ud,ul,ut,fail] = rd.ramses_data(infile,lmax,xc,yc,zc,dx,dy,dz,scalelist[scale])
        
        # Clean exit if the file was not found
        if fail:
            print divider
            return
        
        print "Generating data structure... please wait"
        
        # Create a small information dictionary
        self.info = {"ncells"  : nn      ,\
                     "ncpu"    : ncpu    ,\
                     "ndim"    : ndim    ,\
                     "levelmin": levelmin,\
                     "levelmax": levelmax,\
                     "nstep"   : nstep   ,\
                     "boxsize" : boxsize ,\
                     "time"    : time*ut ,\
                     "ud"      : ud      ,\
                     "ul"      : ul      ,\
                     "ut"      : ut      ,\
                     "center"  : center   }
        
        # This is the master data dictionary. For each entry, the dict has 5 fields.
        # It loops through the list of variables that it got from the file loader.
        # For example, for the density:
        # - data["density"]["values"   ]: holds the values for the density, which have
        #                                 been normalized to cgs using 'ud'
        # - data["density"]["unit"     ]: a string containing the units
        # - data["density"]["label"    ]: a string containing the variable name
        # - data["density"]["operation"]: an operation, if the variable has been
        #                                 derived from other variables. This is empty
        #                                 if it is a core variable.
        # - data["density"]["depth"    ]: a depth which tells you on how many layers
        #                                 of derived variables the current variables
        #                                 depends on.
        self.data = dict()
        list_vars = names.split()
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            self.data[theKey] = dict()
            [norm,uu] = self.get_units(theKey,ud,ul,ut,scale)
            self.data[theKey]["values"   ] = data1[:nn,i]*norm
            self.data[theKey]["unit"     ] = uu
            self.data[theKey]["label"    ] = theKey
            self.data[theKey]["operation"] = ""
            self.data[theKey]["depth"    ] = 0

        # Modifications for coordinates and cell sizes
        self.re_center(center,scale,ndim)
        self.data["dx"]["values"] = self.data["dx"]["values"]/scalelist[scale]
        
        # We now use the 'new_field' function to create commonly used variables
        
        # Magnetic field
        self.new_field(name="B_x",operation="0.5*(B_left_x+B_right_x)",unit="G",label="B_x")
        self.new_field(name="B_y",operation="0.5*(B_left_y+B_right_y)",unit="G",label="B_x")
        self.new_field(name="B_z",operation="0.5*(B_left_z+B_right_z)",unit="G",label="B_x")
        self.new_field(name="B",operation="np.sqrt(B_x**2+B_y**2+B_z**2)",unit="G",label="B")
    
        # Commonly used log quantities
        self.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
        self.new_field(name="log_T",operation="np.log10(temperature)",unit="K",label="log(T)")
        self.new_field(name="log_B",operation="np.log10(B)",unit="g/cm3",label="log(B)")
        
        print infile+" successfully loaded"
        print "The variables are:"
        for key in sorted(self.data):
            print key+" ["+self.data[key]["unit"]+"]"
        print divider
    
    #=======================================================================================
    # The re_center function shifts the coordinates axes around a center. If center="auto"
    # then the function find the cell with the highest density.
    #=======================================================================================
    def re_center(self,center,scale,ndim):
        
        try:
            lc = len(center)
            if center == "auto":
                maxloc = np.argmax(self.data["density"]["values"])
                xc = self.data["x"]["values"][maxloc]
                yc = self.data["y"]["values"][maxloc]
                zc = self.data["z"]["values"][maxloc]
            elif lc == 3:
                xc = center[0]*self.info["boxsize"]
                yc = center[1]*self.info["boxsize"]
                zc = center[2]*self.info["boxsize"]
            else:
                print "Bad center value"
                return
        except TypeError:
            xc = yc = zc = 0.5*self.info["boxsize"]
        self.data["x"]["values"] = (self.data["x"]["values"] - xc)/scalelist[scale]
        self.data["y"]["values"] = (self.data["y"]["values"] - yc)/scalelist[scale]
        if ndim > 2:
            self.data["z"]["values"] = (self.data["z"]["values"] - zc)/scalelist[scale]
            
    #=======================================================================================
    # The new field function is used to create a new data field. Say you want to take the
    # log of the density. You create a new field by calling:
    # mydata.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    # The operation string is then evaluated using the 'eval' function.
    #=======================================================================================
    def new_field(self,name,operation,unit,label):
        
        [op_parsed,depth] = self.parse_operation(operation)
        #print op_parsed
        self.data[name] = dict()
        self.data[name]["values"   ] = eval(op_parsed)
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
        # and after.
        expression = " "+operation+" "
        # Sort the list of variable keys in the order of the longest to the shortest.
        key_list = sorted(self.data.keys(),key=lambda x:len(x),reverse=True)
        #print key_list
        for key in key_list:
            loc = expression.find(key)
            if loc > -1:
                char_before = expression[loc-1]
                char_after  = expression[loc+len(key)]
                bad_before = (char_before.isalpha() or (char_before == "_"))
                bad_after = (char_after.isalpha() or (char_after == "_"))
                #print "|"+char_before+"|  --  |"+char_before+"|"
                if (not bad_before) and (not bad_after):
                    expression = expression.replace(key,"self.data[\""+key+"\"][\"values\"]")
                    max_depth = max(max_depth,self.data[key]["depth"])
            
        return [expression,max_depth]
    
    #=======================================================================================
    #=======================================================================================
    
    def update_values(self,nout=1,lmax=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="cm"):
        
        if nout == -1:
            filelist = sorted(glob.glob("output*"))
            infile = filelist[-1]
        else:
            infile = "output_"+str(nout).zfill(5)
        
        try:
            center += 0
            xc,yc,zc = center[0:3]
        except TypeError:
            xc = yc = zc = 0.5
        
        print divider
        
        [data1,names,nn,ncpu,ndim,levelmin,levelmax,nstep,boxsize,time,ud,ul,ut,fail] = rd.ramses_data(infile,lmax,xc,yc,zc,dx,dy,dz,scalelist[scale])
        
        if fail:
            print divider
            return
        
        print "Updating data structure... please wait"
        
        list_vars = names.split()
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            [norm,uu] = self.get_units(theKey,ud,ul,ut,scale)
            self.data[theKey]["values"   ] = data1[:nn,i]*norm
            self.data[theKey]["unit"     ] = uu
            self.data[theKey]["label"    ] = theKey
            self.data[theKey]["operation"] = ""
        
        # Modifications for coordinates and cell sizes
        self.re_center(center,scale,ndim)
        self.data["dx"]["values"] = self.data["dx"]["values"]/scalelist[scale]
        
        # Now go through the fields and update the values of fields that have a operation attached to them
        # IMPORTANT: this needs to be done in the right order: use the depth key to determine which
        # variables depend on others
        key_list = sorted(self.data.keys(),key=lambda x:self.data[x]["depth"])
        for key in key_list:
            if len(self.data[key]["operation"]) > 0:
                #op_parsed = self.parse_operation(self.data[key]["operation"])
                print "Re-valuating "+key
                self.data[key]["values"] = eval(self.data[key]["operation"])
                
        print "Data successfully updated with values from "+infile
        print divider
        
    #=======================================================================================
    #=======================================================================================

    def plot_histogram(self,var_x,var_y,var_z=None,fname=None,logz=True,axes=None,cmap=None):

        # Parameters
        nx = 101
        ny = 101
        
        datax  = self.data[var_x]["values"]
        datay  = self.data[var_y]["values"]
        xlabel = self.data[var_x]["label"]+" ["+self.data[var_x]["unit"]+"]"
        ylabel = self.data[var_y]["label"]+" ["+self.data[var_y]["unit"]+"]"
                
        xmin = np.amin(datax)
        xmax = np.amax(datax)
        ymin = np.amin(datay)
        ymax = np.amax(datay)
        
        xe = np.linspace(xmin,xmax,nx)
        ye = np.linspace(ymin,ymax,ny)

        z, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe))

        contourz = True
        try:
            dataz = self.data[var_z]["values"]
        except KeyError:
            contourz = False

        if contourz:
            z1, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=dataz)
            z2 = z1/z
            zmin = np.amin(dataz)
            zmax = np.amax(dataz)

        x = np.zeros([nx-1])
        y = np.zeros([ny-1])

        for i in range(nx-1):
            x[i] = 0.5*(xe[i]+xe[i+1])
        for j in range(ny-1):
            y[j] = 0.5*(ye[j]+ye[j+1])

        if logz:
            z = np.log10(z)
        
        # Begin plotting -------------------------------------
        if axes:
            theplot = axes
        else:
            plt.clf()
            plt.subplot(111)
            theplot = plt
            
        cont = theplot.contourf(x,y,z,20,cmap=cmap)
        if contourz:
            over = theplot.contour(x,y,z2,levels=np.arange(zmin,zmax+1),colors='k')
            theplot.clabel(over,inline=1,fmt='%i')
            leg = [over.collections[0]]
            theplot.legend(leg,[self.data[var_z]["label"]],loc=2)
        if axes:
            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)
        else:
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
        if fname:
            plt.savefig(fname,bbox_inches="tight")
        elif axes:
            pass
        else:
            plt.show(block=False)

        return

    #=======================================================================================
    #=======================================================================================

    def plot_slice(self,var="rho",direction="z",vec=None,streamlines=False,fname=None,dx=1.0,dy=0.0,cmap=None,axes=None,resolution=128):
        
        if direction == "z":
            dir_x = "x"
            dir_y = "y"
        elif direction == "y":
            dir_x = "x"
            dir_y = "z"
        elif direction == "x":
            dir_x = "y"
            dir_y = "z"
        else:
            print "Bad direction for slice"
            return
        
        # Make it possible to call with only one size in the arguments
        if dy == 0.0:
            dy = dx

        # Make a guess for slice thickness and iterate if no points are found in the slice
        dz = 0.05*(0.5*(dx+dy))
        no_points = True
        while no_points:
            data1 = self.data[dir_x    ]["values"]
            data2 = self.data[dir_y    ]["values"]
            data3 = self.data[direction]["values"]
            cube  = np.where(np.logical_and(abs(data1) < 0.5*dx,np.logical_and(abs(data2) < 0.5*dy,abs(data3) < 0.5*dz)))
            datax = data1[cube]
            if len(datax) == 0:
                dz = 2.0*dz
            else:
                no_points = False
                
        datay = data2[cube]
        dataz = self.data[var]["values"][cube]
        if vec:
            datau = self.data[vec+"_"+dir_x]["values"][cube]
            datav = self.data[vec+"_"+dir_y]["values"][cube]
        celldx = self.data["dx"]["values"][cube]
        ncells = np.shape(datax)[0]
        
        xmin = -0.5*dx
        xmax =  0.5*dx
        ymin = -0.5*dy
        ymax =  0.5*dy
        
        nx = resolution
        ny = resolution
        dpx = (xmax-xmin)/nx
        dpy = (ymax-ymin)/ny
        
        z1 = np.zeros([ny,nx])
        z2 = np.zeros([ny,nx])
        if vec:
            u1 = np.zeros([ny,nx])
            v1 = np.zeros([ny,nx])
            z3 = np.zeros([ny,nx])
        
        for n in range(ncells):
            x1 = datax[n]-0.5*celldx[n]
            x2 = datax[n]+0.5*celldx[n]
            y1 = datay[n]-0.5*celldx[n]
            y2 = datay[n]+0.5*celldx[n]
            
            ix1 = max(int((x1-xmin)/dpx),0)
            ix2 = min(int((x2-xmin)/dpx),nx-1)
            iy1 = max(int((y1-ymin)/dpy),0)
            iy2 = min(int((y2-ymin)/dpy),ny-1)
                    
            for j in range(iy1,iy2+1):
                for i in range(ix1,ix2+1):
                    z1[j,i] = z1[j,i] + dataz[n]
                    z2[j,i] = z2[j,i] + 1.0
                    if vec:
                        u1[j,i] = u1[j,i] + datau[n]
                        v1[j,i] = v1[j,i] + datav[n]
                        z3[j,i] = z3[j,i] + np.sqrt(datau[n]**2+datav[n]**2)
            
        z = z1/z2
        if vec:
            u = u1/z2
            v = v1/z2
            w = z3/z2
        
        x = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
        y = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,ny)
        iskip = int(0.071*resolution)
        
        xlab = self.data[dir_x]["label"]+" ["+self.data[dir_x]["unit"]+"]"
        ylab = self.data[dir_y]["label"]+" ["+self.data[dir_y]["unit"]+"]"
        zlab = self.data[var  ]["label"]+" ["+self.data[var  ]["unit"]+"]"
        
        # Begin plotting -------------------------------------
        if axes:
            theplot = axes
        else:
            plt.clf()
            plt.subplot(111)
            theplot = plt
        
        cont = theplot.contourf(x,y,z,20,cmap=cmap)
        if axes:
            cbar = plt.colorbar(cont,ax=theplot)
            theplot.set_xlabel(xlab)
            theplot.set_ylabel(ylab)
        else:
            cbar = plt.colorbar(cont)
            plt.xlabel(xlab)
            plt.ylabel(ylab)
        if vec:
            if streamlines:
                if streamlines == "log":
                    w = np.log10(w)
                strm = theplot.streamplot(x,y,u,v,color=w,cmap='Greys')
            else:
                vect = theplot.quiver(x[::iskip],y[::iskip],u[::iskip,::iskip],v[::iskip,::iskip],w[::iskip,::iskip],cmap='Greys',pivot='mid')
        
        cbar.ax.set_ylabel(zlab)
        cbar.ax.yaxis.set_label_coords(-1.0,0.5) 
        
        if fname:
            plt.savefig(fname,bbox_inches="tight")
        elif axes:
            pass
        else:
            plt.show(block=False)
        
        return

    #=======================================================================================
    #=======================================================================================

    def get_units(self,string,ud,ul,ut,scale="cm"):
        if string == "density":
            return [ud,"g/cm3"]
        elif string.startswith("velocity"):
            return [ul/ut,"cm/s"]
        elif string.startswith("B_"):
            return [np.sqrt(4.0*np.pi*ud*(ul/ut)**2),"G"]
        elif string == "thermal_pressure":
            return [ud*((ul/ut)**2),"g/cm/s2"]
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
