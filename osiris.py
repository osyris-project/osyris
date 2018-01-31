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

# Osiris libraries
import config_osiris as conf
import load_ramses_data
import ism_physics
# Python libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata

#=======================================================================================
#=======================================================================================
# TYPES OF CODE OUTPUTS: Define one class per type of snapshot to read
#=======================================================================================
#=======================================================================================

# Ramses data format ===================================================================
class RamsesData(load_ramses_data.LoadRamsesData):
 
    def __init__(self,nout=conf.default_values["nout"],lmax=conf.default_values["lmax"],\
                 center=conf.default_values["center"],dx=conf.default_values["dx"],\
                 dy=conf.default_values["dy"],dz=conf.default_values["dz"],\
                 scale=conf.default_values["scale"],verbose=conf.default_values["verbose"],\
                 path=conf.default_values["path"],variables=conf.default_values["variables"]):
        
        load_ramses_data.LoadRamsesData.__init__(self,nout=nout,lmax=lmax,center=center,\
                 dx=dx,dy=dy,dz=dz,scale=scale,verbose=verbose,path=path,variables=variables)
        
        return
#=======================================================================================

# Heracles data format =================================================================
#=======================================================================================


#=======================================================================================
#=======================================================================================
# OSIRIS PLOTTING FUNCTIONS
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
def plot_histogram(var_x,var_y,scalar=False,image=False,contour=False,fname=None,axes=None,logz=False,\
                   resolution=256,copy=False,title=None,\
                   xmin=None,xmax=None,ymin=None,ymax=None,new_window=False,\
                   update=None,cbar=True,outline=False,scatter=False,summed=False,\
                   clear=True,plot=True,block=False,equal_axes=False,\
                   scalar_args={},image_args={},contour_args={},scatter_args={},outline_args={}):

    # Find parent container of object to plot
    holder = var_x.parent
    
    # Possibility of updating the data from inside the plotting routines
    try:
        update += 0
        holder.update_values(nout=update)
    except TypeError:
        pass

    # Parameters
    nx = resolution+1
    ny = resolution+1
    
    # Get the data values and units
    datax  = var_x.values
    datay  = var_y.values
    xlabel = var_x.label+" ["+var_x.unit+"]"
    ylabel = var_y.label+" ["+var_y.unit+"]"
    default_var = "histo_cell_density"
    
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
    
    # Protect against empty plots if xmin==xmax or ymin==ymax
    if xmin == xmax:
        if xmin == 0.0:
            xmin = -0.1
            xmax = 0.1
        else:
            xmin = xmin - 0.05*abs(xmin)
            xmax = xmax + 0.05*abs(xmax)
    if ymin == ymax:
        if ymin == 0.0:
            ymin = -0.1
            ymax = 0.1
        else:
            ymin = ymin - 0.05*abs(ymin)
            ymax = ymax + 0.05*abs(ymax)
    
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
    
    # Use numpy histogram2d function to make image
    z_scal = z_imag = z_cont = z_outl = False
    empty = True
    if scalar:
        try:
            z1, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=scalar.values)
            if summed:
                z_scal = np.ma.masked_where(z0 == 0.0, z1)
            else:
                with np.errstate(divide="ignore",invalid="ignore"):
                    z_scal = np.ma.masked_where(z0 == 0.0, z1/z0)
        except AttributeError:
            holder.new_field(name=default_var,unit="",label="Number of cells",verbose=True)
            scalar = getattr(holder,default_var)
            z_scal = np.ma.masked_where(z0 == 0.0, z0)
        empty = False
    if image:
        try:
            z1, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=image.values)
            if summed:
                z_imag = np.ma.masked_where(z0 == 0.0, z1)
            else:
                with np.errstate(divide="ignore",invalid="ignore"):
                    z_imag = np.ma.masked_where(z0 == 0.0, z1/z0)
        except AttributeError:
            holder.new_field(name=default_var,unit="",label="Number of cells",verbose=True)
            image = getattr(holder,default_var)
            z_imag = np.ma.masked_where(z0 == 0.0, z0)
        empty = False
    if contour:
        try:
            z1, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=contour.values)
            if summed:
                z_cont = np.ma.masked_where(z0 == 0.0, z1)
            else:
                with np.errstate(divide="ignore",invalid="ignore"):
                    z_cont = np.ma.masked_where(z0 == 0.0, z1/z0)
        except AttributeError:
            holder.new_field(name=default_var,unit="",label="Number of cells",verbose=True)
            contour = getattr(holder,default_var)
            z_cont = np.ma.masked_where(z0 == 0.0, z0)
        empty = False
    if scatter:
        empty = False
    
    # If no variable is requested for z/color dimension, store number of cells by default
    if empty:
        holder.new_field(name=default_var,unit="",label="Number of cells",verbose=True)
        scalar = getattr(holder,default_var)
        z_scal = np.ma.masked_where(z0 == 0.0, z0)
        
    if outline:
        z_outl = z0
        
    if plot:
        render_map(scalar=scalar,image=image,contour=contour,scatter=scatter,x=x,y=y,z_scal=z_scal,    \
                   z_imag=z_imag,z_cont=z_cont,z_outl=z_outl,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fname=fname,          \
                   axes=axes,title=title,new_window=new_window,clear=clear,block=block,          \
                   resolution=resolution,scalar_args=scalar_args,image_args=image_args,    \
                   contour_args=contour_args,scatter_args=scatter_args,equal_axes=equal_axes,x_raw=var_x,y_raw=var_y,\
                   outline=outline,outline_args=outline_args,dir_x=var_x.name,dir_y=var_y.name,sinks=False)
    
    if hasattr(holder,default_var):
        holder.delete_field(default_var)
    
    if copy:
        return x,y,z_scal,z_imag,z_cont,z_outl
    else:
        return

#=======================================================================================
# Plot a 2D slice through the data cube. The arguments are:
# - scalar     : the scalar field to be plotted, e.g. mydata.density
# - image      : the scalar field to be plotted with an image
# - contour    : the scalar field to be plotted with contours
# - vec        : the vector field to be overplotted, e.g. mydata.velocity
# - stream     : the field for streamlines to be overplotted, e.g. mydata.B
# - direction  : the direction normal to the plane of the slice
# - dx         : the x extent of the slice, in units of scale (see data loader)
# - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
# - dz         : the thickness of the slice
# - axes       : if specified, the data is plotted on the specified axes (see demo).
# - resolution : number of pixels in the slice.
# - fname      : if specified, the figure is saved to file.
#=======================================================================================
def plot_slice(scalar=False,image=False,contour=False,vec=False,stream=False,axes=None,\
               direction="z",dx=0.0,dy=0.0,fname=None,title=None,sinks=True,copy=False,\
               origin=[0,0,0],resolution=128,summed=False,new_window=False,update=None,\
               clear=True,plot=True,block=False,interpolation='linear',\
               scalar_args={},image_args={},contour_args={},vec_args={},stream_args={},outline=False,outline_args={}):
    
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
    
    # Possibility of updating the data from inside the plotting routines
    try:
        update += 0
        holder.update_values(nout=update)
    except TypeError:
        pass
    
    # Get slice extent and direction vectors
    dx,dy,box,dir1,dir2,dir3,dir_x,dir_y = get_slice_direction(holder,direction,dx,dy)
    
    # Define equation of a plane
    a_plane = dir1[0]
    b_plane = dir1[1]
    c_plane = dir1[2]
    d_plane = -dir1[0]*origin[0]-dir1[1]*origin[1]-dir1[2]*origin[2]
        
    dist1 = (a_plane*holder.get("x")+b_plane*holder.get("y")+c_plane*holder.get("z")+d_plane) \
          / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)
      
    dist2 = np.sqrt((holder.get("x")-origin[0])**2+(holder.get("y")-origin[1])**2+(holder.get("z")-origin[2])**2) - np.sqrt(3.0)*0.5*holder.get("dx")

    # Select only the cells in contact with the slice., i.e. at a distance less than dx/2
    cube = np.where(np.logical_and(np.abs(dist1) <= 0.5000000001*holder.get("dx"),np.abs(dist2) <= max(dx,dy)*0.5*np.sqrt(2.0)))
    # Project coordinates onto the plane by taking dot product with axes vectors
    coords = np.transpose([holder.get("x")[cube]-origin[0],holder.get("y")[cube]-origin[1],holder.get("z")[cube]-origin[2]])
    datax = np.inner(coords,dir2)
    datay = np.inner(coords,dir3)
    
    # Define slice extent and resolution
    xmin = max(-0.5*dx,box[0])
    xmax = min(xmin+dx,box[1])
    ymin = max(-0.5*dy,box[2])
    ymax = min(ymin+dy,box[3])
    nx   = resolution
    ny   = resolution
    dpx  = (xmax-xmin)/float(nx)
    dpy  = (ymax-ymin)/float(ny)
    x = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
    y = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,ny)
    grid_x, grid_y = np.meshgrid(x, y)
    points = np.transpose([datax,datay])
    
    # Use scipy interpolation function to make image
    z_scal = z_imag = z_cont = u_vect = v_vect = w_vect = u_strm = v_strm = w_strm = 0
    if scalar:
        z_scal = griddata(points,scalar.values[cube] ,(grid_x,grid_y),method=interpolation)
    if image:
        z_imag = griddata(points,image.values[cube]  ,(grid_x,grid_y),method=interpolation)
    if contour:
        z_cont = griddata(points,contour.values[cube],(grid_x,grid_y),method=interpolation)
    if vec:
        if holder.info["ndim"] < 3:
            datau1 = vec.x.values[cube]
            datav1 = vec.y.values[cube]
        else:
            vectors = np.transpose([vec.x.values[cube],vec.y.values[cube],vec.z.values[cube]])
            datau1 = np.inner(vectors,dir2)
            datav1 = np.inner(vectors,dir3)
        u_vect = griddata(points,datau1,(grid_x,grid_y),method=interpolation)
        v_vect = griddata(points,datav1,(grid_x,grid_y),method=interpolation)
        w_vect = griddata(points,np.sqrt(datau1**2+datav1**2),(grid_x,grid_y),method=interpolation)
    if stream:
        if holder.info["ndim"] < 3:
            datau2 = stream.x.values[cube]
            datav2 = stream.y.values[cube]
        else:
            streams = np.transpose([stream.x.values[cube],stream.y.values[cube],stream.z.values[cube]])
            datau2 = np.inner(streams,dir2)
            datav2 = np.inner(streams,dir3)
        u_strm = griddata(points,datau2,(grid_x,grid_y),method=interpolation)
        v_strm = griddata(points,datav2,(grid_x,grid_y),method=interpolation)
        w_strm = griddata(points,np.sqrt(datau2**2+datav2**2),(grid_x,grid_y),method=interpolation)
    
    # Render the map    
    if plot:
        render_map(scalar=scalar,image=image,contour=contour,vec=vec,stream=stream,x=x,y=y,z_scal=z_scal,    \
                   z_imag=z_imag,z_cont=z_cont,u_vect=u_vect,v_vect=v_vect,w_vect=w_vect,u_strm=u_strm,      \
                   v_strm=v_strm,w_strm=w_strm,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fname=fname,          \
                   axes=axes,title=title,sinks=sinks,new_window=new_window,clear=clear,block=block,          \
                   dir_x=dir_x,dir_y=dir_y,resolution=resolution,thePlane=[a_plane,b_plane,c_plane,d_plane], \
                   origin=origin,dir_vecs=[dir1,dir2,dir3],scalar_args=scalar_args,image_args=image_args,    \
                   contour_args=contour_args,vec_args=vec_args,stream_args=stream_args,outline=outline,outline_args=outline_args)
    
    if copy:
        return x,y,z_scal,z_imag,z_cont,u_vect,v_vect,w_vect,u_strm,v_strm,w_strm
    else:
        return

#=======================================================================================
# Plot a column density through the data cube. The arguments are:
# - scalar     : the scalar field to be plotted, e.g. mydata.density
# - image      : the scalar field to be plotted with an image
# - contour    : the scalar field to be plotted with contours
# - dx         : the x extent of the slice, in units of scale (see data loader)
# - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
# - dz         : the thickness of the slice
# - axes       : if specified, the data is plotted on the specified axes (see demo).
# - resolution : number of pixels in the slice.
# - fname      : if specified, the figure is saved to file.
#=======================================================================================
def plot_column_density(scalar=False,image=False,contour=False,vec=False,stream=False,          \
                        direction="z",dx=0.0,dy=0.0,dz=0.0,fname=None,axes=None,title=None,     \
                        origin=[0,0,0],resolution=128,sinks=True,summed=False,copy=False,       \
                        new_window=False,update=None,clear=True,plot=True,block=False,nz=0,     \
                        interpolation="linear",verbose=False,outline=False,outline_args={},\
                        scalar_args={},image_args={},contour_args={},vec_args={},stream_args={}):
        
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
    
    # Possibility of updating the data from inside the plotting routines
    try:
        update += 0
        holder.update_values(nout=update)
    except TypeError:
        pass
        
    # Get direction vectors
    dx,dy,box,dir1,dir2,dir3,dir_x,dir_y = get_slice_direction(holder,direction,dx,dy)
    
    if dz == 0.0:
        dz = max(dx,dy)
    
    # Define equation of a plane
    a_plane = dir1[0]
    b_plane = dir1[1]
    c_plane = dir1[2]
    d_plane = -dir1[0]*origin[0]-dir1[1]*origin[1]-dir1[2]*origin[2]
    
    sqrt3 = np.sqrt(3.0)
    
    # Define slice extent and resolution
    xmin = max(-0.5*dx,box[0])
    xmax = min(xmin+dx,box[1])
    ymin = max(-0.5*dy,box[2])
    ymax = min(ymin+dy,box[3])
    zmin = -0.5*dz
    zmax =  0.5*dz
    nx   = resolution
    ny   = resolution
    if nz == 0:
        nz = resolution
    dpx = (xmax-xmin)/float(nx)
    dpy = (ymax-ymin)/float(ny)
    dpz = (zmax-zmin)/float(nz)
    x = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
    y = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,ny)
    z = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nz)
    grid_x, grid_y = np.meshgrid(x, y)
    
    
    # We now create empty data arrays that will be filled by the cell data
    z_scal = z_imag = z_cont = u_vect = v_vect = w_vect = u_strm = v_strm = w_strm = 0
    if scalar:
        z_scal = np.zeros([ny,nx])
    if image:
        z_imag = np.zeros([ny,nx])
    if contour:
        z_cont = np.zeros([ny,nx])

    iprog = 1
    istep = 10
    
    for iz in range(nz):
        
        # Print progress
        if verbose:
            percentage = int(float(iz)*100.0/float(nz))
            if percentage >= iprog*istep:
                print("%3i%% done" % percentage)
                iprog += 1
    
        dist1 = (a_plane*holder.get("x")+b_plane*holder.get("y")+c_plane*holder.get("z")+d_plane) \
              / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2) - z[iz]
          
        dist2 = np.sqrt((holder.get("x")-origin[0]-z[iz]*dir1[0])**2 + \
                        (holder.get("y")-origin[1]-z[iz]*dir1[1])**2 + \
                        (holder.get("z")-origin[2]-z[iz]*dir1[2])**2) - sqrt3*0.5*holder.get("dx")

        # Select only the cells in contact with the slice., i.e. at a distance less than sqrt(3)*dx/2
        cube = np.where(np.logical_and(np.abs(dist1) <= 0.5000000001*holder.get("dx"),np.abs(dist2) <= max(dx,dy)*0.5*np.sqrt(2.0)))
        #cube = np.where(np.abs(dist1) <= sqrt3*0.5*holder.get("dx"))
        #ncells = np.shape(holder.get("dx")[cube])[0]
        #print ncells
        #celldx = holder.get("dx")[cube]
        # Project coordinates onto the plane by taking dot product with axes vectors
        coords = np.transpose([holder.get("x")[cube]-origin[0]-z[iz]*dir1[0],holder.get("y")[cube]-origin[1]-z[iz]*dir1[1],holder.get("z")[cube]-origin[2]-z[iz]*dir1[2]])
        datax = np.inner(coords,dir2)
        datay = np.inner(coords,dir3)
        points = np.transpose([datax,datay])
        
        if scalar:
            z_scal += griddata(points,scalar.values[cube] ,(grid_x,grid_y),method=interpolation)*dpz*conf.constants[holder.info["scale"]]
        if image:
            z_imag += griddata(points,image.values[cube]  ,(grid_x,grid_y),method=interpolation)*dpz*conf.constants[holder.info["scale"]]
        if contour:
            z_cont += griddata(points,contour.values[cube],(grid_x,grid_y),method=interpolation)*dpz*conf.constants[holder.info["scale"]]
    
    # Render the map    
    if plot:
        render_map(scalar=scalar,image=image,contour=contour,vec=vec,stream=stream,x=x,y=y,z_scal=z_scal,    \
                   z_imag=z_imag,z_cont=z_cont,u_vect=u_vect,v_vect=v_vect,w_vect=w_vect,u_strm=u_strm,      \
                   v_strm=v_strm,w_strm=w_strm,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fname=fname,          \
                   axes=axes,title=title,sinks=sinks,new_window=new_window,clear=clear,block=block,          \
                   dir_x=dir_x,dir_y=dir_y,resolution=resolution,thePlane=[a_plane,b_plane,c_plane,d_plane], \
                   origin=origin,dir_vecs=[dir1,dir2,dir3],scalar_args=scalar_args,image_args=image_args,    \
                   contour_args=contour_args,vec_args=vec_args,stream_args=stream_args,outline=outline,outline_args=outline_args)
    
    if copy:
        return x,y,z_scal,z_imag,z_cont,u_vect,v_vect,w_vect,u_strm,v_strm,w_strm
    else:
        return

#=======================================================================================
# Compute a vector perpendicular to the input vector
#=======================================================================================
def perpendicular_vector(v):

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0:
        raise ValueError("zero-vector")
    
    if v[2] == 0:
        return [-v[1],v[0],0]
    else:
        return [1.0, 1.0, -1.0 * (v[0] + v[1]) / v[2]]

#=======================================================================================
# Separate arguments to plotting functions between Osiris specific arguments and
# Matplotlib arguments.
#=======================================================================================
def parse_arguments(args,args_osiris,args_plot):
    
    # Save a copy so we can exclude these parameters specific to osiris from the ones
    # to be sent to the matplotlib routines.
    args_osiris["cb_format"] = None
    ignore = set(args_osiris.keys())
    # Now we go through the arguments taken from the function call - scalar_args - and
    # add them to the osiris arguments.
    for key in args.keys():
        args_osiris[key] = args[key]
    # Define colorbar scaling
    cmap = args_osiris["cmap"]
    norm = None
    # Default contour levels
    try:
        clevels = np.linspace(args_osiris["vmin"],args_osiris["vmax"],args_osiris["nc"])
    except TypeError:
        clevels = [0.0,1.0]
    # Perform log normalization if it is required
    try:
        if cmap.startswith("log") or cmap.endswith("log"):
            #cmap_save = cmap
            cmap = cmap.replace("log","")
            chars = [",",":",";"," "]
            for ch in chars:
                cmap = cmap.replace(ch,"")
            if len(cmap) == 0:
                cmap = conf.default_values["colormap"] # cmap_save
            #norm = LogNorm()
            norm = LogNorm(vmin=args_osiris["vmin"],vmax=args_osiris["vmax"])
            clevels = np.logspace(np.log10(args_osiris["vmin"]),np.log10(args_osiris["vmax"]),args_osiris["nc"])
            args_osiris["cb_format"] = "%.1e"
    except AttributeError:
        pass
    
    # We now define default parameters for the matplotlib function
    keylist = args_plot.keys()
    if "levels" in keylist:
        args_plot["levels"]=clevels
    if "cmap" in keylist:
        args_plot["cmap"]=cmap
    if "norm" in keylist:
        args_plot["norm"]=norm
    # Then run through the vec_args, adding them to the plotting arguments, but
    # ignoring all osiris specific arguments
    keys = set(args.keys())
    for key in keys.difference(ignore):
        args_plot[key] = args[key]
    
    return

#=======================================================================================
# Find direction vectors for slice
#=======================================================================================
def get_slice_direction(holder,direction,dx,dy):
    
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
        vel    = np.vstack((holder.get("velocity_x")[sphere],holder.get("velocity_y")[sphere],holder.get("velocity_z")[sphere])).T
        #vel    = holder.get("velocity")[sphere]
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
    
    box = [boxmin_x,boxmax_x,boxmin_y,boxmax_y]

    return dx,dy,box,dir1,dir2,dir3,dir_x,dir_y

#=======================================================================================
# Use matplotlib to plot histogram, slice or column density maps
#=======================================================================================
def render_map(scalar=False,image=False,contour=False,scatter=False,vec=False,stream=False,outline=False,x=0,y=0,  \
               z_scal=0,z_imag=0,z_cont=0,z_outl=0,u_vect=0,v_vect=0,w_vect=0,u_strm=0,v_strm=0,\
               w_strm=0,fname=None,axes=None,title=None,sinks=True,new_window=False,   \
               clear=True,block=False,xmin=0,xmax=0,ymin=0,ymax=0,dir_x="x",dir_y="y", \
               resolution=128,scalar_args={},image_args={},contour_args={},vec_args={},\
               stream_args={},scatter_args={},outline_args={},dz=0,thePlane=0,origin=[0,0,0],\
               dir_vecs=0,x_raw=None,y_raw=None,equal_axes=True):
    
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
    elif scatter:
        holder = x_raw.parent
    else:
        print("Nothing to render.")
        return
    
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
    
    if scalar:
        
        # Round off AMR levels to integers
        if scalar.label == "level" or scalar.label == "cpu":
            z_scal = np.around(z_scal)
        # Parse scalar plot arguments
        scalar_args_osiris = {"vmin":np.nanmin(z_scal),"vmax":np.nanmax(z_scal),"cbar":True,"cbax":None,"cmap":conf.default_values["colormap"],"nc":21}
        scalar_args_plot = {"levels":1,"cmap":1,"norm":1}
        parse_arguments(scalar_args,scalar_args_osiris,scalar_args_plot)
        contf = theAxes.contourf(x,y,z_scal,**scalar_args_plot)
        if scalar_args_osiris["cbar"]:
            scb = plt.colorbar(contf,ax=theAxes,cax=scalar_args_osiris["cbax"],format=scalar_args_osiris["cb_format"])
            scb.ax.set_ylabel(scalar.label+(" ["+scalar.unit+"]" if len(scalar.unit) > 0 else ""))
            scb.ax.yaxis.set_label_coords(-1.1,0.5)
            
    if image:
        
        # Round off AMR levels to integers
        if image.label == "level" or image.label == "cpu":
            z_imag = np.around(z_imag)
        # Here we define a set of default parameters
        image_args_osiris = {"vmin":np.nanmin(z_imag),"vmax":np.nanmax(z_imag),"cbar":True,"cbax":None,"cmap":conf.default_values["colormap"],"nc":21}
        # cmap and norm are just dummy arguments to tell the parsing function that they are required
        image_args_plot = {"cmap":1,"norm":1,"interpolation":"none","origin":"lower"}
        parse_arguments(image_args,image_args_osiris,image_args_plot)
        img = theAxes.imshow(z_imag,extent=[xmin,xmax,ymin,ymax],**image_args_plot)
        if image_args_osiris["cbar"]:
            icb = plt.colorbar(img,ax=theAxes,cax=image_args_osiris["cbax"],format=image_args_osiris["cb_format"])
            icb.ax.set_ylabel(image.label+(" ["+image.unit+"]" if len(image.unit) > 0 else ""))
            icb.ax.yaxis.set_label_coords(-1.1,0.5)
            
    if contour:
        
        # Round off AMR levels to integers
        if contour.label == "level" or contour.label == "cpu":
            z_cont = np.around(z_cont)
        # Here we define a set of default parameters
        contour_args_osiris = {"vmin":np.nanmin(z_cont),"vmax":np.nanmax(z_cont),"cbar":False,"cbax":None,\
                               "cmap":conf.default_values["colormap"],"nc":21,"label":False,"fmt":"%1.3f"}
        # levels, cmap and norm are just dummy arguments to tell the parsing function that they are required
        contour_args_plot = {"levels":1,"cmap":1,"norm":1,"zorder":10,"linestyles":"solid"}
        parse_arguments(contour_args,contour_args_osiris,contour_args_plot)
        cont = theAxes.contour(x,y,z_cont,**contour_args_plot)
        if contour_args_osiris["label"]:
            theAxes.clabel(cont,inline=1,fmt=contour_args_osiris["fmt"])
        if contour_args_osiris["cbar"]:
            ccb = plt.colorbar(cont,ax=theAxes,cax=contour_args_osiris["cbax"],format=contour_args_osiris["cb_format"])
            ccb.ax.set_ylabel(contour.label+(" ["+contour.unit+"]" if len(contour.unit) > 0 else ""))
            ccb.ax.yaxis.set_label_coords(-1.1,0.5)
    
    if scatter:
        try:
            vmin = np.nanmin(scatter.values)
            vmax = np.nanmax(scatter.values)
            scbar = True
            scmap = conf.default_values["colormap"]
        except AttributeError:
            vmin = vmax = None
            scbar = False
            scmap = None
        scatter_args_osiris = {"iskip":1,"cmap":scmap,"vmin":vmin,"vmax":vmax,"cbar":scbar,"cbax":None,"nc":21}
        scatter_args_plot = {"cmap":1,"marker":".","c":"b","edgecolor":"None","s":20,"norm":1}
        parse_arguments(scatter_args,scatter_args_osiris,scatter_args_plot)
        # Check if a variable is given as a color
        try:
            scatter_args_plot["c"] = scatter.values[::scatter_args_osiris["iskip"]]
        except AttributeError:
            pass
        scat = theAxes.scatter(x_raw.values[::scatter_args_osiris["iskip"]],y_raw.values[::scatter_args_osiris["iskip"]],**scatter_args_plot)
        if scatter_args_osiris["cbar"]:
            rcb = plt.colorbar(scat,ax=theAxes,cax=scatter_args_osiris["cbax"],format=scatter_args_osiris["cb_format"])
            rcb.ax.set_ylabel(scatter.label+(" ["+scatter.unit+"]" if len(scatter.unit) > 0 else ""))
            rcb.ax.yaxis.set_label_coords(-1.1,0.5)
    
    if outline:        
        outline_args_plot = {"levels":[np.nanmin(z_outl)],"colors":"grey"}
        for key in outline_args.keys():
            outline_args_plot[key] = outline_args[key]
        outl = theAxes.contour(x,y,z_outl, **outline_args_plot)
    
    # Plot vector field
    if vec:
        
        vec_args_osiris = {"vskip":int(0.047*resolution),"vscale":np.nanmax(w_vect),"vsize":15.0,"vkey":True,"vkey_pos":[0.70,-0.08],\
                           "cbar":False,"cbax":None,"vmin":np.nanmin(w_vect),"vmax":np.nanmax(w_vect),"nc":21,"cmap":None}
        vec_args_plot = {"cmap":1,"pivot":"mid","color":"w","norm":None}
        parse_arguments(vec_args,vec_args_osiris,vec_args_plot)
        vskip = vec_args_osiris["vskip"]
        if not "scale" in vec_args_plot.keys():
            vec_args_plot["scale"] = vec_args_osiris["vsize"]*vec_args_osiris["vscale"]
        if vec_args_plot["cmap"]:
            vect = theAxes.quiver(x[::vskip],y[::vskip],u_vect[::vskip,::vskip],v_vect[::vskip,::vskip],\
                                  w_vect[::vskip,::vskip],**vec_args_plot)
            if vec_args_osiris["cbar"]:
                vcb = plt.colorbar(vect,ax=theAxes,cax=vec_args_osiris["cbax"],orientation="horizontal",format=vec_args_osiris["cb_format"])
                vcb.ax.set_xlabel(vec.label+(" ["+vec.unit+"]" if len(vec.unit) > 0 else ""))
        else:
            vect = theAxes.quiver(x[::vskip],y[::vskip],u_vect[::vskip,::vskip],v_vect[::vskip,::vskip],\
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
        stream_args_osiris = {"cbar":False,"cbax":None,"sskip":1,"vmin":np.nanmin(w_strm),"vmax":np.nanmax(w_strm),"nc":21,"cmap":None}
        stream_args_plot = {"cmap":1,"color":"w","norm":None}
        parse_arguments(stream_args,stream_args_osiris,stream_args_plot)
        sskip = stream_args_osiris["sskip"]
        if stream_args_plot["cmap"]:
            stream_args_plot["color"]=w_strm[::sskip,::sskip]
        strm = theAxes.streamplot(x[::sskip],y[::sskip],u_strm[::sskip,::sskip],v_strm[::sskip,::sskip],**stream_args_plot)
        if stream_args_osiris["cbar"]:
            scb = plt.colorbar(strm.lines,ax=theAxes,cax=stream_args_osiris["cbax"],orientation="horizontal",format=stream_args_osiris["cb_format"])
            scb.ax.set_xlabel(stream.label+(" ["+stream.unit+"]" if len(stream.unit) > 0 else ""))
        
    
    if holder.info["nsinks"] > 0 and sinks:
        dx = xmax-xmin
        if dz == 0:
            thickness = 0.05*dx
        else:
            thickness = 0.5*dz
        dist = (thePlane[0]*holder.sinks["x"]+thePlane[1]*holder.sinks["y"]+thePlane[2]*holder.sinks["z"]+thePlane[3]) \
               / np.sqrt(thePlane[0]**2 + thePlane[1]**2 + thePlane[2]**2)
        sinkcoords = np.transpose([holder.sinks["x"]-origin[0],holder.sinks["y"]-origin[1],holder.sinks["z"]-origin[2]])
        sink_x = np.inner(sinkcoords,dir_vecs[1])
        sink_y = np.inner(sinkcoords,dir_vecs[2])
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
    
    #if clear:
        #theAxes.set_xlim([xmin,xmax])
        #theAxes.set_ylim([ymin,ymax])
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
    
    if equal_axes:
        theAxes.set_aspect("equal")

    if fname:
        plt.savefig(fname,bbox_inches="tight")
    elif axes:
        pass
    else:
        plt.show(block=block)
        
    return

#=======================================================================================
# Interpolate data at any given point in the whole 3D domain
#=======================================================================================
def interpolate(field,points):
    
    holder = field.parent
    
    try:
        hashTable = holder.hash_table
    except AttributeError:
        print("A hash table is needed to perform interpolations")
        holder.create_hash_table()
    
    points[:,0] = ((points[:,0] + holder.info["xc"])/holder.info["boxsize_scaled"])
    points[:,1] = ((points[:,1] + holder.info["yc"])/holder.info["boxsize_scaled"])
    points[:,2] = ((points[:,2] + holder.info["zc"])/holder.info["boxsize_scaled"])

    npoints = np.shape(points)[0]    
    ilevl = holder.info["levelmax"]
    values = np.zeros([npoints])
    for ip in range(npoints):
        not_found = True
        loop_count = 0
        while not_found:
            l = max(min(ilevl+((-1)**loop_count)*int((loop_count+1)/2),holder.info["levelmax"]),0)
            loop_count += 1
            dxcell = 0.5**l
            igrid = int(points[ip,0]/dxcell)
            jgrid = int(points[ip,1]/dxcell)
            kgrid = int(points[ip,2]/dxcell)
            theHash = str(igrid)+','+str(jgrid)+','+str(kgrid)+','+str(l)
            try:
                icell = holder.hash_table[theHash]
                ilevl = l
                not_found = False
            except KeyError:
                pass
        
        cube = dict()
        dmax = 0.0
        cube[theHash] = dict()
        cube[theHash]["vars"] = field.values[holder.hash_table[theHash]]
        cube[theHash]["dist"] = np.sqrt((points[ip,0]-((igrid+0.5)*dxcell))**2 + \
                                        (points[ip,1]-((jgrid+0.5)*dxcell))**2 + \
                                        (points[ip,2]-((kgrid+0.5)*dxcell))**2)
        dmax = max(dmax,cube[theHash]["dist"])
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if i == j == k == 1:
                        pass
                    else:
                        ii = igrid-1+i
                        jj = jgrid-1+j
                        kk = kgrid-1+k
                        theHash = str(ii)+','+str(jj)+','+str(kk)+','+str(ilevl)
                        try:
                            neighbour = holder.hash_table[theHash]
                            cube[theHash] = dict()
                            cube[theHash]["vars"] = field.values[holder.hash_table[theHash]]
                            cube[theHash]["dist"] = np.sqrt((points[ip,0]-((ii+0.5)*dxcell))**2 + \
                                                            (points[ip,1]-((jj+0.5)*dxcell))**2 + \
                                                            (points[ip,2]-((kk+0.5)*dxcell))**2)
                            dmax = max(dmax,cube[theHash]["dist"])
                        except KeyError:
                            theHash = str(2*ii)+','+str(2*jj)+','+str(2*kk)+','+str(ilevl+1)
                            try:
                                neighbour = holder.hash_table[theHash]
                                for i1 in range(2):
                                    for j1 in range(2):
                                        for k1 in range(2):
                                            theHash = str(2*ii+i1)+','+str(2*jj+j1)+','+str(2*kk+k1)+','+str(ilevl+1)
                                            cube[theHash] = dict()
                                            cube[theHash]["vars"] = field.values[holder.hash_table[theHash]]
                                            cube[theHash]["dist"] = np.sqrt((points[ip,0]-((2*ii+i1+0.5)*dxcell*0.5))**2 + \
                                                                            (points[ip,1]-((2*jj+j1+0.5)*dxcell*0.5))**2 + \
                                                                            (points[ip,2]-((2*kk+k1+0.5)*dxcell*0.5))**2)
                                            dmax = max(dmax,cube[theHash]["dist"])
                            except KeyError:
                                theHash = str(int(float(ii)/2.0))+','+str(int(float(jj)/2.0))+','+str(int(float(kk)/2.0))+','+str(ilevl-1)
                                try:
                                    neighbour = holder.hash_table[theHash]
                                    cube[theHash] = dict()
                                    cube[theHash]["vars"] = field.values[holder.hash_table[theHash]]
                                    cube[theHash]["dist"] = np.sqrt((points[ip,0]-((int(float(ii)/2.0)+0.5)*dxcell*2.0))**2 + \
                                                                    (points[ip,1]-((int(float(jj)/2.0)+0.5)*dxcell*2.0))**2 + \
                                                                    (points[ip,2]-((int(float(kk)/2.0)+0.5)*dxcell*2.0))**2)
                                    dmax = max(dmax,cube[theHash]["dist"])
                                except KeyError:
                                    print("Neighbour not found",igrid,jgrid,kgrid,i,j,k)
                        
        # Compute inverse distance weighting
        result  = 0.0
        weights = 0.0
        for key in cube.keys():
            #w = (0.1-1.0)/dmax * cube[key]["dist"] + 1.0
            w = 1.0 / (np.exp(15.0*(cube[key]["dist"]/dmax-0.5)) + 1.0)
            weights += w
            result  += w*cube[key]["vars"]
        
        values[ip] = result/weights

    return values

#=======================================================================================
# Write RAMSES data to VTK file
#=======================================================================================
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
