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
    if var_x.kind == "vector":
        datax = np.linalg.norm(var_x.values,axis=1)
    else:
        datax = var_x.values
    if var_y.kind == "vector":
        datay = np.linalg.norm(var_y.values,axis=1)
    else:
        datay  = var_y.values
    xlabel = var_x.label+" ["+var_x.unit+"]"
    ylabel = var_y.label+" ["+var_y.unit+"]"
    if var_z:
        if var_z.kind == "vector":
            dataz = np.linalg.norm(var_z.values,axis=1)
        else:
            dataz  = var_z.values
            zlabel = var_z.label
    
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
            with np.errstate(divide="ignore"):
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
                        
        theAxes.set_xlabel(xlabel)
        theAxes.set_ylabel(ylabel)
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

#=======================================================================================
# Plot a 2D slice through the data cube. The arguments are:
# - var        : the key for the variable to be plotted, e.g. "density" or "log_rho"
# - direction  : the direction normal to the plane of the slice
# - vec        : the vector field to be overplotted. For velocity, one should supply
#                "velocity" as input and the routine will search for "velocity_x" and
#                "velocity_y" in the variable fields.
# - stream     : the field for streamlines to be overplotted. For B field, one should
#                supply "B" as input and the routine will search for "B_x" and
#                "B_y" in the variable fields.
# - fname      : if specified, the figure is saved to file.
# - dx         : the x extent of the slice, in units of scale (see data loader)
# - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
# - cmap       : the colormap
# - axes       : if specified, the data is plotted on the specified axes (see demo).
# - resolution : number of pixels in the slice.
#=======================================================================================
def plot_slice(scalar=False,vec=False,stream=False,direction="z",fname=None,\
               dx=0.0,dy=0.0,dz=0.0,cmap=conf.default_values["colormap"],axes=None,\
               nc=20,new_window=False,sinks=True,update=None,vmin=None,vmax=None,\
               title=None,cbar=True,cbax=None,clear=True,plot=True,block=False,\
               origin=[0,0,0],summed=False,image=False,resolution=128,copy=False,\
               contour=False,scalar_args={},image_args={},contour_args={},\
               vec_args={},stream_args={}):
    
    ## Possibility of updating the data from inside the plotting routines
    #try:
        #update += 0
        #self.update_values(nout=update)
    #except TypeError:
        #pass
    
    # Find parent container of object to plot
    if scalar:
        holder = scalar.parent
    elif image:
        holder = image.parent
    elif contour:
        holder = contour.parent
    elif vec:
        holder = vec.parent
    elif stream:
        holder = stream.parent
    else:
        print("Nothing to plot.")
        return
        
    if holder.info["ndim"] < 2:
        print("Cannot plot slice from 1D data. Exiting...")
        return
    
    # List of directions
    dir_list = {"x" : ["y","z"], "y" : ["x","z"], "z" : ["x","y"], "auto" : ["x","y"], "auto:top" : ["x","y"], "auto:side" : ["x","z"]}
    
    
    
    
    # Set dx to whole box if not specified
    boxmin_x = np.nanmin(holder.get(dir_list.get(direction,["x","y"])[0]))
    boxmax_x = np.nanmax(holder.get(dir_list.get(direction,["x","y"])[0]))
    boxmin_y = np.nanmin(holder.get(dir_list.get(direction,["x","y"])[1]))
    boxmax_y = np.nanmax(holder.get(dir_list.get(direction,["x","y"])[1]))
    if dx+dy == 0.0:
        dx = boxmax_x - boxmin_x
        dy = boxmax_y - boxmin_y
    elif dx == 0.0:
        dx = dy
    # Make it possible to call with only one size in the arguments
    if dy == 0.0:
        dy = dx
    
    # Define x,y directions depending on the input direction
    if direction[0]=="[" and direction[-1]=="]":
        dir_x = "x"
        dir_y = "y"
        dir1 = eval(direction)
        dir2 = perpendicular_vector(dir1)
        dir3 = np.cross(dir1,dir2)
    elif direction.startswith("auto"):
        params = direction.split(":")
        if len(params) == 1:
            view = "top"
        else:
            view = params[1]
        if len(params) < 3:
            sphere_rad = 0.5*((np.nanmax(holder.get("x"))-np.nanmin(holder.get("x"))) if dx == 0.0 else dx)
        else:
            sphere_rad = float(params[2])
        dir_x = "x"
        dir_y = "y"
        # Compute angular momentum vector
        sphere = np.where(holder.get("r") < sphere_rad)
        pos    = np.vstack((holder.get("x")[sphere],holder.get("y")[sphere],holder.get("z")[sphere])*holder.get("mass")[sphere]).T
        #vel    = np.vstack((holder.get("velocity_x")[sphere],holder.get("velocity_y")[sphere],holder.get("velocity_z")[sphere])).T
        vel    = holder.get("velocity")[sphere]
        AngMom = np.sum(np.cross(pos,vel),axis=0)
        if view == "top":
            dir1 = AngMom
            dir2 = perpendicular_vector(dir1) # [1.0, 1.0, -1.0 * (dir1[0] + dir1[1]) / dir1[2]]
            dir3 = np.cross(dir1,dir2)
        elif view == "side":
            # Choose a vector perpendicular to the angular momentum vector
            dir3 = AngMom
            dir1 = perpendicular_vector(dir3) # [1.0, 1.0, -1.0 * (dir3[0] + dir3[1]) / dir3[2]]
            dir2 = np.cross(dir1,dir3)
        norm1 = np.linalg.norm(dir1)
        print("Normal slice vector: [%.5e,%.5e,%.5e]" % (dir1[0]/norm1,dir1[1]/norm1,dir1[2]/norm1))
    elif ((direction == "x") or (direction == "y") or (direction == "z")):
        [dir_x,dir_y] = dir_list[direction]
        dir1 = [int(direction=="x"),int(direction=="y"),int(direction=="z")]
        dir2 = [int(direction=="y" or direction=="z"),int(direction=="x"),0]
        dir3 = [0,int(direction=="z"),int(direction=="x" or direction=="y")]
    elif holder.info["ndim"]==2:
        dir1 = [0,0,1]
        dir2 = [1,0,0]
        dir3 = [0,1,0]
        dir_x = "x"
        dir_y = "y"
    else:
        print("Bad direction for slice")
        return
    
    norm1 = np.linalg.norm(dir1)
    norm2 = np.linalg.norm(dir2)
    norm3 = np.linalg.norm(dir3)
    dir1 = dir1 / norm1
    dir2 = dir2 / norm2
    dir3 = dir3 / norm3
    
    # Define equation of a plane
    a_plane = dir1[0]
    b_plane = dir1[1]
    c_plane = dir1[2]
    d_plane = -dir1[0]*origin[0]-dir1[1]*origin[1]-dir1[2]*origin[2]
    
    sqrt3 = np.sqrt(3.0)
    
    dist = (a_plane*holder.get("x")+b_plane*holder.get("y")+c_plane*holder.get("z")+d_plane) \
         / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)

    # Select only the cells in contact with the slice., i.e. at a distance less than sqrt(3)*dx/2
    cube = np.where(abs(dist) <= sqrt3*0.5*holder.get("dx")+0.5*dz)
    ncells = np.shape(holder.get("dx")[cube])[0]
    celldx = holder.get("dx")[cube]
    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = np.transpose([holder.get("x")[cube]-origin[0],holder.get("y")[cube]-origin[1],holder.get("z")[cube]-origin[2]])
    datax = np.inner(coords,dir2)
    datay = np.inner(coords,dir3)
    
    if scalar:
        if scalar.kind == "vector":
            dataz1 = np.linalg.norm(scalar.values[cube,:],axis=1)
        else:
            dataz1 = scalar.values[cube]
    
    if image:
        if image.kind == "vector":
            dataz2 = np.linalg.norm(image.values[cube,:],axis=1)
        else:
            dataz2 = image.values[cube]
    
    if contour:
        if contour.kind == "vector":
            dataz3 = np.linalg.norm(contour.values[cube,:],axis=1)
        else:
            dataz3 = contour.values[cube]
    
    # Now project vectors and streamlines using the same method
    if vec:
        if vec.kind == "scalar":
            print("Warning: cannot make vectors out of scalar field.")
        else:
            if holder.info["ndim"] < 3:
                datau1 = vec.values[cube][0]
                datav1 = vec.values[cube][1]
            else:
                #print vec.values[cube]
                #vectors = np.transpose([holder.get(vec+"_x")[cube],holder.get(vec+"_y")[cube],vz])
                #datau1 = np.inner(vec.values[cube,:],dir2)
                #datav1 = np.inner(vec.values[cube,:],dir3)
                datau1 = np.inner(vec.values[cube],dir2)
                datav1 = np.inner(vec.values[cube],dir3)
        #print np.shape(datau1)
        #print np.shape(vec.values[cube,:]),np.shape(dir2)
    if stream:
        if stream.kind == "scalar":
            print("Warning: cannot make streamlines out of scalar field.")
        else:
            if holder.info["ndim"] < 3:
                datau2 = stream.values[cube][0]
                datav2 = stream.values[cube][1]
            else:
                #vectors = np.transpose([holder.get(vec+"_x")[cube],holder.get(vec+"_y")[cube],vz])
                datau2 = np.inner(stream.values[cube],dir2)
                datav2 = np.inner(stream.values[cube],dir3)
        
        #try:
            #sz = holder.get(stream+"_z")[cube]
        #except KeyError:
            #sz = holder.get(stream+"_x")[cube]*0.0
        #streams = np.transpose([holder.get(stream+"_x")[cube],holder.get(stream+"_y")[cube],sz])
        #datau2 = np.inner(streams,dir2)
        #datav2 = np.inner(streams,dir3)
    
    # Define slice extent and resolution
    xmin = max(-0.5*dx,boxmin_x)
    xmax = min(xmin+dx,boxmax_x)
    ymin = max(-0.5*dy,boxmin_y)
    ymax = min(ymin+dy,boxmax_y)
    nx   = resolution
    ny   = resolution
    dpx  = (xmax-xmin)/float(nx)
    dpy  = (ymax-ymin)/float(ny)
    
    # We now create empty data arrays that will be filled by the cell data
    za = np.zeros([ny,nx])
    zb = np.zeros([ny,nx])
    zc = np.zeros([ny,nx])
    zd = np.zeros([ny,nx])
    if vec:
        u1 = np.zeros([ny,nx])
        v1 = np.zeros([ny,nx])
        z1 = np.zeros([ny,nx])
    if stream:
        u2 = np.zeros([ny,nx])
        v2 = np.zeros([ny,nx])
        z2 = np.zeros([ny,nx])
    
    # Loop through all data cells and find extent covered by the current cell size
    for n in range(ncells):
        x1 = datax[n]-0.5*celldx[n]*sqrt3
        x2 = datax[n]+0.5*celldx[n]*sqrt3
        y1 = datay[n]-0.5*celldx[n]*sqrt3
        y2 = datay[n]+0.5*celldx[n]*sqrt3
        
        # Find the indices of the slice pixels which are covered by the current cell
        ix1 = max(int((x1-xmin)/dpx),0)
        ix2 = min(int((x2-xmin)/dpx),nx-1)
        iy1 = max(int((y1-ymin)/dpy),0)
        iy2 = min(int((y2-ymin)/dpy),ny-1)
        
        # Fill in the slice pixels with data
        for j in range(iy1,iy2+1):
            for i in range(ix1,ix2+1):
                za[j,i] = za[j,i] + celldx[n]
                if scalar:
                    zb[j,i] = zb[j,i] + dataz1[n]*celldx[n]
                if image:
                    zc[j,i] = zc[j,i] + dataz2[n]*celldx[n]
                if contour:
                    zd[j,i] = zd[j,i] + dataz3[n]*celldx[n]
                
                if vec:
                    u1[j,i] = u1[j,i] + datau1[n]*celldx[n]
                    v1[j,i] = v1[j,i] + datav1[n]*celldx[n]
                    z1[j,i] = z1[j,i] + np.sqrt(datau1[n]**2+datav1[n]**2)*celldx[n]
                if stream:
                    u2[j,i] = u2[j,i] + datau2[n]*celldx[n]
                    v2[j,i] = v2[j,i] + datav2[n]*celldx[n]
                    z2[j,i] = z2[j,i] + np.sqrt(datau2[n]**2+datav2[n]**2)*celldx[n]
    
    # Compute z averages
    if summed:
        if scalar:
            z_scal = np.ma.masked_where(za == 0.0, zb)
        if image:
            z_imag = np.ma.masked_where(za == 0.0, zc)
        if contour:
            z_cont = np.ma.masked_where(za == 0.0, zd)
        
        if vec:
            u1 = np.ma.masked_where(za == 0.0, u1)
            v1 = np.ma.masked_where(za == 0.0, v1)
            w1 = np.ma.masked_where(za == 0.0, z1)
        if stream:
            u2 = np.ma.masked_where(za == 0.0, u2)
            v2 = np.ma.masked_where(za == 0.0, v2)
            w2 = np.ma.masked_where(za == 0.0, z2)
    else:
        if scalar:
            z_scal = np.ma.masked_where(za == 0.0, zb/za)
        if image:
            z_imag = np.ma.masked_where(za == 0.0, zc/za)
        if contour:
            z_cont = np.ma.masked_where(za == 0.0, zd/za)
        if vec:
            u1 = np.ma.masked_where(za == 0.0, u1/za)
            v1 = np.ma.masked_where(za == 0.0, v1/za)
            w1 = np.ma.masked_where(za == 0.0, z1/za)
        if stream:
            u2 = np.ma.masked_where(za == 0.0, u2/za)
            v2 = np.ma.masked_where(za == 0.0, v2/za)
            w2 = np.ma.masked_where(za == 0.0, z2/za)
        
    # Define cell centers for filled contours
    x = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
    y = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,ny)
    
    
    
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
        
        # If vmin and vmax are specified as a normal argument and not in a dict,
        # there should only be one of scalar, image or contour
        score = 0
        if scalar:
            score += 1
        if image:
            score += 1
        if contour:
            score += 1

        #if score > 1:
            #print("vmin and vmax limits cannot be applied to more than one scalar field.")
            #print("Please use scalar_args={}, image_args={}, contour_args={}.")
            #return
            
        
        if scalar:
            
            # Round off AMR levels to integers
            if scalar.label == "level":
                z_scal = np.around(z_scal)
            
            scalar_args_osiris  = {"vmin"   : int(0.047*resolution),
                                "vscale"  : np.nanmax(w1),
                                "vsize"   : 15.0,
                                "vkey"    : True,
                                "vkey_pos": [0.70,-0.08],
                                "cbar"    : False,
                                "cbax"    : None}
            # Save a copy so we can exclude these parameters specific to osiris from the ones
            # to be sent to the matplotlib quiver routine.
            ignore = set(vec_args_osiris.keys())
            # Now we go through the arguments taken from the function call - vec_args - and
            # add them to the osiris arguments.
            for key in vec_args.keys():
                vec_args_osiris[key] = vec_args[key]
            # We now define default parameters for the quiver function
            vec_args_plot = {"cmap":None,"pivot":"mid","scale":vec_args_osiris["vsize"]*vec_args_osiris["vscale"],"color":"w","norm":None}
            # Then run through the vec_args, adding them to the plotting arguments, but
            # ignoring all osiris specific arguments
            keys = set(vec_args.keys())
            for key in keys.difference(ignore):
                vec_args_plot[key] = vec_args[key]
            
            # Define colorbar limits
            try:
                zmin += 0
            except TypeError:
                zmin = np.nanmin(z_scal)
            try:
                zmax += 0
            except TypeError:
                zmax = np.nanmax(z_scal)
            
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
            
            slice_args_plot = {"levels":clevels,"cmap":cmap,"norm":norm}
            for key in slice_args.keys():
                slice_args_plot[key] = slice_args[key]
            contf = theAxes.contourf(x,y,z_scal,**slice_args_plot)
            if cbar:
                cb = plt.colorbar(contf,ax=theAxes,cax=cbax,format=cb_format)
                z_scal_lab = scalar.label
                if len(scalar.unit) > 0:
                    z_scal_lab += " ["+scalar.unit+"]"
                cb.ax.set_ylabel(z_scal_lab)
                cb.ax.yaxis.set_label_coords(-1.1,0.5)
            
        if image:
            # Define colorbar limits
            try:
                zmin += 0
            except TypeError:
                zmin = np.nanmin(z_imag)
            try:
                zmax += 0
            except TypeError:
                zmax = np.nanmax(z_imag)
            
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
            
            image_args_plot = {"interpolation":"none","origin":"lower","cmap":cmap,"norm":norm}
            for key in image_args.keys():
                image_args_plot[key] = image_args[key]
            img = theAxes.imshow(z_imag,extent=[xmin,xmax,ymin,ymax],vmin=zmin,vmax=zmax,**image_args_plot)
            
            if cbar:
                cb = plt.colorbar(img,ax=theAxes,cax=cbax,format=cb_format)
                z_imag_lab = image.label
                if len(image.unit) > 0:
                    z_imag_lab += " ["+image.unit+"]"
                cb.ax.set_ylabel(z_imag_lab)
                cb.ax.yaxis.set_label_coords(-1.1,0.5)
                
        if contour:
            # Define colorbar limits
            try:
                zmin += 0
            except TypeError:
                zmin = np.nanmin(z_cont)
            try:
                zmax += 0
            except TypeError:
                zmax = np.nanmax(z_cont)
            
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
            
            contour_args_plot = {"levels":clevels,"cmap":cmap,"norm":norm,"zorder":10}
            clabel_args = {"label":False,"fmt":"%1.3f"}
            ignore = set(clabel_args.keys())
            keys = set(contour_args.keys())
            for key in keys.difference(ignore):
                contour_args_plot[key] = contour_args[key]
            for key in contour_args.keys():
                clabel_args[key] = contour_args[key]
            cont = theAxes.contour(x,y,z_cont,**contour_args_plot)
            if clabel_args["label"]:
                theAxes.clabel(cont,inline=1,fmt=clabel_args["fmt"])
            if cbar:
                cb = plt.colorbar(cont,ax=theAxes,cax=cbax,format=cb_format)
                z_cont_lab = contour.label
                if len(contour.unit) > 0:
                    z_cont_lab += " ["+contour.unit+"]"
                cb.ax.set_ylabel(z_cont_lab)
                cb.ax.yaxis.set_label_coords(-1.1,0.5)
        #theAxes.set_xlabel(xlab)
        #theAxes.set_ylabel(ylab)
        
        # Plot vector field
        if vec:
            # Here we define a set of default parameters
            vec_args_osiris  = {"vskip"   : int(0.047*resolution),
                                "vscale"  : np.nanmax(w1),
                                "vsize"   : 15.0,
                                "vkey"    : True,
                                "vkey_pos": [0.70,-0.08],
                                "cbar"    : False,
                                "cbax"    : None}
            # Save a copy so we can exclude these parameters specific to osiris from the ones
            # to be sent to the matplotlib quiver routine.
            ignore = set(vec_args_osiris.keys())
            # Now we go through the arguments taken from the function call - vec_args - and
            # add them to the osiris arguments.
            for key in vec_args.keys():
                vec_args_osiris[key] = vec_args[key]
            # We now define default parameters for the quiver function
            vec_args_plot = {"cmap":None,"pivot":"mid","scale":vec_args_osiris["vsize"]*vec_args_osiris["vscale"],"color":"w","norm":None}
            # Then run through the vec_args, adding them to the plotting arguments, but
            # ignoring all osiris specific arguments
            keys = set(vec_args.keys())
            for key in keys.difference(ignore):
                vec_args_plot[key] = vec_args[key]
            # We are now ready to plot the vectors. Note that we need two different calls if
            # a colormap is used for the vectors.
            vskip = vec_args_osiris["vskip"]
            vcmap = vec_args_plot["cmap"]
            if vcmap:
                if vcmap.startswith("log") or vcmap.endswith("log"):
                    vcmap = vcmap.replace("log","")
                    chars = [",",":",";"," "]
                    for ch in chars:
                        vcmap = vcmap.replace(ch,"")
                    if len(vcmap) == 0:
                        vcmap = conf.default_values["colormap"]
                    vec_args_plot["cmap"] = vcmap
                    vec_args_plot["norm"] = LogNorm(vmin=np.nanmin(w1[::vskip,::vskip]),vmax=np.nanmax(w1[::vskip,::vskip]))
                    vcb_format = "%.1e"
                else:
                    vcb_format = None
                
                vect = theAxes.quiver(x[::vskip],y[::vskip],u1[::vskip,::vskip],v1[::vskip,::vskip],\
                                      w1[::vskip,::vskip],**vec_args_plot)
                if vec_args_osiris["cbar"]:
                    vcb = plt.colorbar(vect,ax=theAxes,cax=vec_args_osiris["cbax"],orientation="horizontal",format=vcb_format)
                    vcb.ax.set_xlabel(vec.label+" ["+vec.unit+"]")
            else:
                vect = theAxes.quiver(x[::vskip],y[::vskip],u1[::vskip,::vskip],v1[::vskip,::vskip],\
                                      **vec_args_plot)
            
            # Plot the scale of the vectors under the axes
            unit_u = vec.unit
            if vec_args_osiris["vkey"]:
                theAxes.quiverkey(vect,vec_args_osiris["vkey_pos"][0],vec_args_osiris["vkey_pos"][1],\
                                  vec_args_osiris["vscale"],"%.2f [%s]" % (vec_args_osiris["vscale"],\
                                  unit_u),labelpos="E",labelcolor="k",coordinates="axes", color="k", \
                                  zorder=100)

        if stream:
            # Here we define a set of default parameters
            stream_args_osiris = {"cbar":False,"cbax":None,"sskip":1}
            # Save a copy so we can exclude these parameters specific to osiris from the ones
            # to be sent to the matplotlib quiver routine.
            ignore = set(stream_args_osiris.keys())
            # Now we go through the arguments taken from the function call - stream_args - and
            # add them to the osiris arguments.
            for key in stream_args.keys():
                stream_args_osiris[key] = stream_args[key]
            # We now define default parameters for the streamplot function
            stream_args_plot = {"cmap":None,"color":"w","norm":None}
            # Then run through the stream_args, adding them to the plotting arguments, but
            # ignoring all osiris specific arguments
            keys = set(stream_args.keys())
            for key in keys.difference(ignore):
                stream_args_plot[key] = stream_args[key]
            # We are now ready to plot the streamlines.
            sskip = stream_args_osiris["sskip"]
            scmap = stream_args_plot["cmap"]
            if scmap:
                if scmap.startswith("log") or scmap.endswith("log"):
                    scmap = scmap.replace("log","")
                    chars = [",",":",";"," "]
                    for ch in chars:
                        scmap = scmap.replace(ch,"")
                    if len(scmap) == 0:
                        scmap = conf.default_values["colormap"]
                    stream_args_plot["cmap"] = scmap
                    stream_args_plot["norm"] = LogNorm(vmin=np.nanmin(w2[::sskip,::sskip]),vmax=np.nanmax(w2[::sskip,::sskip]))
                    scb_format = "%.1e"
                else:
                    scb_format = None
                stream_args_plot["color"]=w2[::sskip,::sskip]
            strm = theAxes.streamplot(x[::sskip],y[::sskip],u2[::sskip,::sskip],v2[::sskip,::sskip],**stream_args_plot)
            if stream_args_osiris["cbar"]:
                    scb = plt.colorbar(strm.lines,ax=theAxes,cax=stream_args_osiris["cbax"],orientation="horizontal",format=scb_format)
                    scb.ax.set_xlabel(stream.label+" ["+stream.unit+"]")
        
        if holder.info["nsinks"] > 0 and sinks:
            if dz == 0.0:
                thickness = 0.05*dx
            else:
                thickness = 0.5*dz
            dist = (a_plane*holder.sinks["x"]+b_plane*holder.sinks["y"]+c_plane*holder.sinks["z"]+d_plane) / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)
            sinkcoords = np.transpose([holder.sinks["x"]-origin[0],holder.sinks["y"]-origin[1],holder.sinks["z"]-origin[2]])
            sink_x = np.inner(sinkcoords,dir2)
            sink_y = np.inner(sinkcoords,dir3)
            subset = np.where(np.logical_and(dist <= thickness,np.logical_and(np.absolute(sink_x) <= 0.5*dx,np.absolute(sink_y) <= 0.5*dx)))
            srad = np.maximum(holder.sinks["radius"][subset],np.full(len(subset),dx*0.01))
            xy = np.array([sink_x[subset],sink_y[subset]]).T
            patches = [plt.Circle(cent, size) for cent, size in zip(xy, srad)]
            coll = matplotlib.collections.PatchCollection(patches, facecolors='w',edgecolors="k",linewidths=2,alpha=0.7)
            theAxes.add_collection(coll)
            
        try:
            title += ""
            theAxes.set_title(title)
        except TypeError:
            theAxes.set_title("Time = %.3f %s" % (holder.info["time"]/conf.constants[conf.default_values["time_unit"]],conf.default_values["time_unit"]))
        
        if clear:
            theAxes.set_xlim([xmin,xmax])
            theAxes.set_ylim([ymin,ymax])
        #else:
            #theAxes.set_xlim([min(theAxes.get_xlim()[0],xmin),max(theAxes.get_xlim()[1],xmax)])
            #theAxes.set_ylim([min(theAxes.get_ylim()[0],ymin),max(theAxes.get_ylim()[1],ymax)])
        
        # Define axes labels
        xlab = getattr(holder,dir_x).label
        if len(getattr(holder,dir_x).unit) > 0:
            xlab += " ["+getattr(holder,dir_x).unit+"]"
        ylab = getattr(holder,dir_y).label
        if len(getattr(holder,dir_y).unit) > 0:
            ylab += " ["+getattr(holder,dir_y).unit+"]"
        theAxes.set_xlabel(xlab)
        theAxes.set_ylabel(ylab)
        
        theAxes.set_aspect("equal")

        if fname:
            plt.savefig(fname,bbox_inches="tight")
        elif axes:
            pass
        else:
            plt.show(block=block)
    
    if copy:
        if vec and stream:
            return x,y,z,u1,v1,w1,u2,v2,w2
        elif vec:
            return x,y,z,u1,v1,w1
        elif stream:
            return x,y,z,u2,v2,w2
        else:
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



def perpendicular_vector(v):

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0:
        raise ValueError("zero-vector")
    
    if v[2] == 0:
        return [-v[1],v[0],0]
    else:
        return [1.0, 1.0, -1.0 * (v[0] + v[1]) / v[2]]
