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

import config_osiris as conf
import load_ramses_data
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib.colors import LogNorm

#=======================================================================================
#np.seterr(divide="ignore",invalid="ignore") # Ignore divide by zero warnings
#=======================================================================================

# Define one class per type of snapshot to read ============================================

# Ramses data format =======================================================================
class RamsesData(load_ramses_data.LoadRamsesData):
 
    def __init__(self,nout=conf.default_values["nout"],lmax=conf.default_values["lmax"],\
                 center=conf.default_values["center"],dx=conf.default_values["dx"],\
                 dy=conf.default_values["dy"],dz=conf.default_values["dz"],\
                 scale=conf.default_values["scale"],verbose=conf.default_values["verbose"],\
                 path=conf.default_values["path"],variables=conf.default_values["variables"]):
        
        load_ramses_data.LoadRamsesData.__init__(self,nout=nout,lmax=lmax,center=center,\
                 dx=dx,dy=dy,dz=dz,scale=scale,verbose=verbose,path=path,variables=variables)
        
        return

# Heracles data format =====================================================================

#=======================================================================================
#=======================================================================================
#=======================================================================================
#=======================================================================================



#=======================================================================================
# Plot a 2D histogram with two variables as input. This is used for instance to plot the
# temperature as a function of density for every cell in the mesh. The input arguments
# are:
# - var_x: a string containing the key for the variable along the x axis, e.g. "log_rho"
# - var_y: a string containing the key for the variable along the y axis, e.g. "log_T"
# - var_z: a string containing the key for a 3rd variable whose contours as overlayed
# - fname: if specified, the figure is saved to file
# - logz : if True, the colormap is logarithmic
# - axes : if specified, the data is plotted on the specified axes (see demo).
# - cmap : the colormap
# - resolution: the data is binned in a 2D matrix of size 'resolution' 
#=======================================================================================
def plot_histogram(var_x,var_y,var_z=None,contour=False,fname=None,axes=None,\
                   cmap=conf.default_values["colormap"],resolution=256,copy=False,\
                   xmin=None,xmax=None,ymin=None,ymax=None,nc=20,new_window=False,\
                   update=None,cbar=True,outline=False,scatter=False,summed=False,\
                   clear=True,plot=True,block=False,zmin=None,zmax=None,cbax=None,\
                   histogram_args={},scatter_args={},contour_args={},outline_args={}):

    ## Possibility of updating the data from inside the plotting routines
    #try:
        #update += 0
        #self.update_values(nout=update)
    #except TypeError:
        #pass

    # Parameters
    nx = resolution+1
    ny = resolution+1
    
    # Get the data values and units
    datax  = var_x.values
    datay  = var_y.values
    #xlabel = self.data[var_x]["label"]+" ["+self.data[var_x]["unit"]+"]"
    #ylabel = self.data[var_y]["label"]+" ["+self.data[var_y]["unit"]+"]"
    if var_z:
        dataz  = var_z.values
        #zlabel = self.data[var_z]["label"]
    
    # Define plotting range
    autoxmin = False
    autoxmax = False
    autoymin = False
    autoymax = False
    try:
        xmin += 0
    except TypeError:
        datax[np.isneginf(datax)] = np.inf
        xmin = np.nanmin(datax)
        autoxmin = True
    try:
        xmax += 0
    except TypeError:
        datax[np.isinf(datax)] = np.NINF
        xmax = np.nanmax(datax)
        autoxmax = True
    try:
        ymin += 0
    except TypeError:
        datay[np.isneginf(datay)] = np.inf
        ymin = np.nanmin(datay)
        autoymin = True
    try:
        ymax += 0
    except TypeError:
        datay[np.isinf(datay)] = np.NINF
        ymax = np.nanmax(datay)
        autoymax = True
    
    dx = xmax-xmin
    dy = ymax-ymin
    if autoxmin:
        xmin = xmin - 0.05*dx
    if autoxmax:
        xmax = xmax + 0.05*dx
    if autoymin:
        ymin = ymin - 0.05*dy
    if autoymax:
        ymax = ymax + 0.05*dy
    
    # Construct some edge specifiers for the histogram2d function call
    xe = np.linspace(xmin,xmax,nx)
    ye = np.linspace(ymin,ymax,ny)
    # Call the numpy histogram2d function
    z0, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe))
    # In the contour plots, x and y are the centers of the cells, instead of the edges.
    x = np.zeros([nx-1])
    y = np.zeros([ny-1])
    for i in range(nx-1):
        x[i] = 0.5*(xe[i]+xe[i+1])
    for j in range(ny-1):
        y[j] = 0.5*(ye[j]+ye[j+1])
    
    if var_z:
        z1, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=dataz)
        if summed:
            z = np.ma.masked_where(z0 == 0.0, z1)
        else:
            z = np.ma.masked_where(z0 == 0.0, z1/z0)
        #zlabel = self.data[var_z]["label"]+" ["+self.data[var_z]["unit"]+"]"
    else:
        z = np.ma.masked_where(z0 == 0.0, z0)
        zlabel = "Number of cells"
    
    # Begin plotting -------------------------------------
    if plot:
        if axes:
            theAxes = axes
        elif new_window:
            plt.figure()
            plt.subplot(111)
            theAxes = plt.gca()
        else:
            if clear:
                plt.clf()
            plt.subplot(111)
            theAxes = plt.gca()
        
        # Define contour limits
        try:
            zmin += 0
        except TypeError:
            zmin = np.nanmin(z)
        try:
            zmax += 0
        except TypeError:
            zmax = np.nanmax(z)
        
        if cmap.startswith("log") or cmap.endswith("log"):
            cmap = cmap.replace("log","")
            chars = [",",":",";"," "]
            for ch in chars:
                cmap = cmap.replace(ch,"")
            if len(cmap) == 0:
                cmap = conf.default_values["colormap"]
            norm = LogNorm()
            clevels = np.logspace(np.log10(zmin),np.log10(zmax),nc)
            cb_format = "%.1e"
        else:
            norm = None
            clevels = np.linspace(zmin,zmax,nc)
            cb_format = None
        
        if scatter:
            scatter_args_osiris = {"iskip":1}
            # Save a copy so we can exclude these parameters specific to osiris from the ones
            # to be sent to the matplotlib quiver routine.
            ignore = set(scatter_args_osiris.keys())
            # Now we go through the arguments taken from the function call - scatter_args - and
            # add them to the osiris arguments.
            for key in scatter_args.keys():
                scatter_args_osiris[key] = scatter_args[key]
            # We now define default parameters for the scatter function
            scatter_args_plot = {"cmap":cmap,"marker":".","c":"b","edgecolor":"None","s":20,"norm":norm}
            # Then run through the vec_args, adding them to the plotting arguments, but
            # ignoring all osiris specific arguments
            keys = set(scatter_args.keys())
            for key in keys.difference(ignore):
                scatter_args_plot[key] = scatter_args[key]
            iskip = scatter_args_osiris["iskip"]
            # Now we can start plotting
            if var_z:
                scatter_args_plot["c"] = dataz[::iskip]
            cont = theAxes.scatter(datax[::iskip],datay[::iskip],**scatter_args_plot)
        elif contour:
            # Run through arguments
            contour_args_plot = {"levels":clevels,"cmap":cmap,"norm":norm}
            clabel_args = {"label":False,"fmt":"%1.3f"}
            ignore = set(clabel_args.keys())
            keys = set(contour_args.keys())
            for key in keys.difference(ignore):
                contour_args_plot[key] = contour_args[key]
            for key in contour_args.keys():
                clabel_args[key] = contour_args[key]
            cont = theAxes.contour(x,y,z,nc,**contour_args_plot)
            if clabel_args["label"]:
                theAxes.clabel(cont,inline=1,fmt=clabel_args["fmt"])
        else:
            # Run through arguments
            histogram_args_plot = {"levels":clevels,"cmap":cmap,"norm":norm}
            for key in histogram_args.keys():
                histogram_args_plot[key] = histogram_args[key]
            # First plot the filled colour contours
            cont = theAxes.contourf(x,y,z,nc,**histogram_args_plot)
                    
        if outline:
            outline_args_plot = {"levels":[1.0],"colors":"grey"}
            for key in outline_args.keys():
                outline_args_plot[key] = outline_args[key]
            outl = theAxes.contour(x,y,z0,**outline_args_plot)
        
        if ((var_z or (not scatter)) and (cbar)):
            cb = plt.colorbar(cont,ax=theAxes,cax=cbax,format=cb_format)
            cb.ax.set_ylabel(zlabel)
            cb.ax.yaxis.set_label_coords(-1.1,0.5)
                        
        #theAxes.set_xlabel(xlabel)
        #theAxes.set_ylabel(ylabel)
        if clear:
            theAxes.set_xlim([xmin,xmax])
            theAxes.set_ylim([ymin,ymax])
        else:
            theAxes.set_xlim([min(theAxes.get_xlim()[0],xmin),max(theAxes.get_xlim()[1],xmax)])
            theAxes.set_ylim([min(theAxes.get_ylim()[0],ymin),max(theAxes.get_ylim()[1],ymax)])
        if fname:
            plt.savefig(fname,bbox_inches="tight")
        elif axes:
            pass
        else:
            plt.show(block=block)

    if copy:
        return x,y,z
    else:
        return















def to_vtk(holder,fname="osiris_data.vtu",variables=False):
    
    try:
        from scipy.spatial import Delaunay
    except ImportError:
        print("Scipy Delaunay library not found. This is needed for VTK output. Exiting.")

    # Print status
    if not fname.endswith(".vtu"):
        fname += ".vtu"
    print("Writing data to VTK file: "+fname)
    
    # Coordinates ot RAMSES cell centers
    #print getattr(getattr( holder,'x'),'values')
    points = np.array([holder.get("x"),holder.get("y"),holder.get("z")]).T
    
    # Compute Delaunay tetrahedralization from cell nodes
    # Note that this step can take a lot of time!
    ncells = holder.info["ncells"]
    print("Computing Delaunay mesh with %i points." % ncells)
    print("This may take some time...")
    tri = Delaunay(points)
    ntetra = np.shape(tri.simplices)[0]
    nverts = ntetra*4
    print("Delaunay mesh with %i tetrahedra complete." % ntetra)

    # Create list of variables by grouping x,y,z components together
    nvarmax = len(holder.data.keys())
    n_components = [] #np.zeros([nvarmax],dtype=np.int32)
    varlist = []
    varnames = []
    for key in holder.data.keys():
        if key.endswith("_x") or key.endswith("_y") or key.endswith("_z"):
            rawkey = key[:-2]
            try:
                k = len(holder.get(rawkey+"_x"))+len(holder.get(rawkey+"_y"))+len(holder.get(rawkey+"_z"))
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
    #for key in holder.data.keys():
        if n_components[i] == 3:
            celldata = np.ravel(np.array([holder.get(varlist[i][0]),holder.get(varlist[i][1]),holder.get(varlist[i][2])]).T)
        else:
            celldata = holder.get(varlist[i][0])
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