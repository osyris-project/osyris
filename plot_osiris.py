import numpy as np
import osiris_common
import config_osiris as conf
import matplotlib.pyplot as plt

#=======================================================================================
# This is a dummy class which gives access to the plotting functions to the other
# classes through inheritance.
#=======================================================================================
class OsirisData(osiris_common.OsirisCommon):
     
    def __init__(self):
                
        return
                
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
    def plot_histogram(self,var_x,var_y,var_z=None,var_c=None,fname=None,logz=False,axes=None,\
                       cmap=conf.default_values["colormap"],resolution=256,copy=False,xmin=None,\
                       xmax=None,ymin=None,ymax=None,nc=20,new_window=False,evol=False,update=None,\
                       outline=False,scatter=False,marker=".",iskip=1,color="b",summed=False,\
                       cbar=True,clear=True,plot=True,block=False):

        # Possibility of updating the data from inside the plotting routines
        try:
            update += 0
            self.update_values(nout=update)
        except TypeError:
            pass

        # Parameters
        nx = resolution+1
        ny = resolution+1
        
        # Get the data values and units
        datax  = self.get(var_x)
        datay  = self.get(var_y)
        xlabel = self.data[var_x]["label"]+" ["+self.data[var_x]["unit"]+"]"
        ylabel = self.data[var_y]["label"]+" ["+self.data[var_y]["unit"]+"]"
        if var_z:
            dataz  = self.get(var_z)
            zlabel = self.data[var_z]["label"]
        
        # Define plotting range
        autoxmin = False
        autoxmax = False
        autoymin = False
        autoymax = False
        
        try:
            xmin += 0
        except TypeError:
            xmin = np.amin(datax)
            autoxmin = True
        try:
            xmax += 0
        except TypeError:
            xmax = np.amax(datax)
            autoxmax = True
        try:
            ymin += 0
        except TypeError:
            ymin = np.amin(datay)
            autoymin = True
        try:
            ymax += 0
        except TypeError:
            ymax = np.amax(datay)
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
        
        if (outline or (not scatter)):
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
        
        if scatter:
            xs = datax[::iskip]
            ys = datay[::iskip]
            if var_z:
                zs = dataz[::iskip]
            else:
                zs = color
        else:
        
            if var_z:
                z1, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=dataz)
                if summed:
                    z = np.ma.masked_where(z0 == 0.0, z1)
                else:
                    z = np.ma.masked_where(z0 == 0.0, z1/z0)
                if logz:
                    z = np.log10(z)
                    zlabel = "log("+zlabel+")"
                zlabel += " ["+self.data[var_z]["unit"]+"]"
            else:
                z = np.ma.masked_where(z0 == 0.0, np.log10(z0))
                zlabel = "log(Number of cells)"
                
            if var_c:
                datac = self.get(var_c)
                z2, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe),weights=datac)
                c = np.ma.masked_where(z0 == 0.0, z2/z0)
        
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
            
            if scatter:
                cont = theAxes.scatter(xs,ys,c=zs,marker=marker,edgecolor='None',cmap=cmap)
            else:
                # First plot the filled colour contours
                cont = theAxes.contourf(x,y,z,nc,cmap=cmap)
            
                # If var_c is specified, overlay black contours
                if var_c:
                    over = theAxes.contour(x,y,c,colors="k")
                    theAxes.clabel(over,inline=1)
                    leg = [over.collections[0]]
                    theAxes.legend(leg,[self.data[var_c]["label"]],loc=2)
            
            if outline:
                outl = theAxes.contour(x,y,z0,levels=[1.0],colors="grey")
            
            if ((var_z or (not scatter)) and (cbar)):
                cb = plt.colorbar(cont,ax=theAxes)
                cb.ax.set_ylabel(zlabel)
                cb.ax.yaxis.set_label_coords(-1.1,0.5)
            
            # Plot evolution (this is quite specific to star formation)
            if evol:
                f = open(evol,"a")
                imax = np.argmax(datax)
                f.write("%.14e  %.14e  %.14e\n" % (self.info["time"],datax[imax],datay[imax]))
                f.close()
                datat = np.loadtxt(evol)
                if np.shape(np.shape(datat))[0] > 1:
                    theAxes.plot(datat[:,1],datat[:,2],color="k",lw=3)
                
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
    # - streamlines: if true, streamlines are drawn for the vector fields instead of arrows.
    #                In addition, if you set streamlines="log", the coloring of the
    #                streamlines will be logarithmic.
    # - fname      : if specified, the figure is saved to file.
    # - dx         : the x extent of the slice, in units of scale (see data loader)
    # - dy         : the y extent of the slice, in units of scale. If not specified, dy = dx
    # - cmap       : the colormap
    # - axes       : if specified, the data is plotted on the specified axes (see demo).
    # - resolution : number of pixels in the slice.
    #=======================================================================================
    def plot_slice(self,var="density",direction="z",vec=False,stream=False,fname=None,\
                   dx=0.0,dy=0.0,dz=0.0,cmap=conf.default_values["colormap"],axes=None,\
                   nc=20,new_window=False,vcmap=False,scmap=False,sinks=True,update=None,\
                   zmin=None,zmax=None,extend="neither",vscale=None,vsize=15.0,title=None,\
                   vcolor="w",scolor="w",vkey_pos=[0.70,-0.08],cbar=True,cbax=None,clear=True,
                   vskip=None,vkey=True,plot=True,center=False,block=False,origin=[0,0,0],
                   summed=False,vsum=False,ssum=False,logz=False,image=False,resolution=128,\
                   copy=False):
        
        # Possibility of updating the data from inside the plotting routines
        try:
            update += 0
            self.update_values(nout=update)
        except TypeError:
            pass
        
        # List of directions
        dir_list = {"x" : ["y","z"], "y" : ["x","z"], "z" : ["x","y"]}
        
        # Set dx to whole box if not specified
        boxmin_x = np.amin(self.get(dir_list.get(direction)[0]))
        boxmax_x = np.amax(self.get(dir_list.get(direction)[0]))
        boxmin_y = np.amin(self.get(dir_list.get(direction)[1]))
        boxmax_y = np.amax(self.get(dir_list.get(direction)[1]))
        if dx+dy == 0.0:
            dx = boxmax_x - boxmin_x
            dy = boxmax_y - boxmin_y
        elif dx == 0.0:
            dx = dy
        # Make it possible to call with only one size in the arguments
        if dy == 0.0:
            dy = dx
        
        # Define x,y directions depending on the input direction
        if len(direction) == 3:
            dir_x = "x"
            dir_y = "y"
            dir1 = direction
        elif direction.startswith("auto"):
            params = direction.split(":")
            if len(params) == 1:
                view = "top"
            else:
                view = params[1]
            if len(params) < 3:
                sphere_rad = 0.5*((np.amax(self.get("x"))-np.amin(self.get("x"))) if dx == 0.0 else dx)
            else:
                sphere_rad = float(params[2])
            dir_x = "x"
            dir_y = "y"
            # Compute angular momentum vector
            sphere = np.where(self.get("r") < sphere_rad)
            pos    = np.vstack((self.get("x")[sphere],self.get("y")[sphere],self.get("z")[sphere])*self.get("mass")[sphere]).T
            vel    = np.vstack((self.get("velocity_x")[sphere],self.get("velocity_y")[sphere],self.get("velocity_z")[sphere])).T
            AngMom = np.sum(np.cross(pos,vel),axis=0)
            if view == "top":
                dir1 = AngMom
            elif view == "side":
                # Choose a vector perpendicular to the angular momentum vector
                dir1 = osiris_common.perpendicular_vector(AngMom)
            print("Normal slice vector: [%.5e,%.5e,%.5e]" % (dir1[0],dir1[1],dir1[2]))
        elif ((direction == "x") or (direction == "y") or (direction == "z")):
            [dir_x,dir_y] = dir_list[direction]
            dir1 = [int(direction=="x"),int(direction=="y"),int(direction=="z")]
        elif self.info["ndim"]==2:
            dir1 = [0,0,1]
            dir_x = "x"
            dir_y = "y"
        else:
            print("Bad direction for slice")
            return
        
        # Choose 2 vectors normal to the direction n and normal to each other
        dir2 = osiris_common.perpendicular_vector(dir1)
        dir3 = np.cross(dir1,dir2)
        
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
        
        dist = (a_plane*self.get("x")+b_plane*self.get("y")+c_plane*self.get("z")+d_plane) \
             / np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)

        # Select only the cells in contact with the slice., i.e. at a distance less than sqrt(3)*dx/2
        cube = np.where(abs(dist) <= sqrt3*0.5*self.get("dx")+0.5*dz)
                        
        dataz = self.get(var)[cube]
        ncells = np.shape(dataz)[0]
        celldx = self.get("dx")[cube]
        
        # Project coordinates onto the plane by taking dot product with axes vectors
        coords = np.transpose([self.get("x")[cube]-origin[0],self.get("y")[cube]-origin[1],self.get("z")[cube]-origin[2]])
        datax = np.inner(coords,dir2)
        datay = np.inner(coords,dir3)
        # Now project vectors and streamlines using the same method
        if vec:
            try:
                vz = self.get(vec+"_z")[cube]
            except KeyError:
                vz = self.get(vec+"_x")[cube]*0.0
            vectors = np.transpose([self.get(vec+"_x")[cube],self.get(vec+"_y")[cube],vz])
            datau1 = np.inner(vectors,dir2)
            datav1 = np.inner(vectors,dir3)
        if stream:
            try:
                sz = self.get(stream+"_z")[cube]
            except KeyError:
                sz = self.get(stream+"_x")[cube]*0.0
            streams = np.transpose([self.get(stream+"_x")[cube],self.get(stream+"_y")[cube],sz])
            datau2 = np.inner(streams,dir2)
            datav2 = np.inner(streams,dir3)
        
        # Define slice extent and resolution
        xmin = max(-0.5*dx,boxmin_x)
        xmax = xmin + dx
        ymin = max(-0.5*dy,boxmin_y)
        ymax = ymin + dy
        nx   = resolution
        ny   = resolution
        dpx  = (xmax-xmin)/nx
        dpy  = (ymax-ymin)/ny
        
        # We now create empty data arrays that will be filled by the cell data
        za = np.zeros([ny,nx])
        zb = np.zeros([ny,nx])
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
                    za[j,i] = za[j,i] + dataz[n]*celldx[n]
                    zb[j,i] = zb[j,i] + celldx[n]
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
            z = np.ma.masked_where(zb == 0.0, za)
        else:
            z = np.ma.masked_where(zb == 0.0, za/zb)
        if logz:
            z = np.log10(z)
        if vec:
            if vsum:
                u1 = np.ma.masked_where(zb == 0.0, u1)
                v1 = np.ma.masked_where(zb == 0.0, v1)
                w1 = np.ma.masked_where(zb == 0.0, z1)
            else:
                u1 = np.ma.masked_where(zb == 0.0, u1/zb)
                v1 = np.ma.masked_where(zb == 0.0, v1/zb)
                w1 = np.ma.masked_where(zb == 0.0, z1/zb)
        if stream:
            if ssum:
                u2 = np.ma.masked_where(zb == 0.0, u2)
                v2 = np.ma.masked_where(zb == 0.0, v2)
                w2 = np.ma.masked_where(zb == 0.0, z2)
            else:
                u2 = np.ma.masked_where(zb == 0.0, u2/zb)
                v2 = np.ma.masked_where(zb == 0.0, v2/zb)
                w2 = np.ma.masked_where(zb == 0.0, z2/zb)
        
        # Define cell centers for filled contours
        x = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
        y = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,ny)
        
        # Define axes labels
        xlab = self.data[dir_x]["label"]
        if len(self.data[dir_x]["unit"]) > 0:
            xlab += " ["+self.data[dir_x]["unit"]+"]"
        ylab = self.data[dir_y]["label"]
        if len(self.data[dir_y]["unit"]) > 0:
            ylab += " ["+self.data[dir_y]["unit"]+"]"
        zlab = self.data[var  ]["label"]
        if len(self.data[var  ]["unit"]) > 0:
            zlab += " ["+self.data[var  ]["unit"]+"]"
        
        # Define colorbar limits
        need_levels = False
        try:
            zmin += 0
            need_levels = True
        except TypeError:
            zmin = np.amin(z)
        try:
            zmax += 0
            need_levels = True
        except TypeError:
            zmax = np.amax(z)
        if need_levels:
            clevels = np.linspace(zmin,zmax,nc)
        else:
            clevels = None
        
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
            
            if image:
                cont = theAxes.imshow(z,extent=[xmin,xmax,ymin,ymax],cmap=cmap,vmin=zmin,vmax=zmax,interpolation="none",origin="lower")
            else:
                cont = theAxes.contourf(x,y,z,nc,levels=clevels,cmap=cmap,extend=extend)
            if cbar:
               cb = plt.colorbar(cont,ax=theAxes,cax=cbax)
               cb.ax.set_ylabel(zlab)
               cb.ax.yaxis.set_label_coords(-1.1,0.5)
            theAxes.set_xlabel(xlab)
            theAxes.set_ylabel(ylab)
            
            if vec:
                try:
                    vskip += 0
                except TypeError:
                    vskip = int(0.071*resolution)
                
                try:
                    vscale += 0
                except TypeError:
                    vscale = np.amax(w1[::vskip,::vskip])
                
                if vcmap:
                    vect = theAxes.quiver(x[::vskip],y[::vskip],u1[::vskip,::vskip],v1[::vskip,::vskip],\
                                          w1[::vskip,::vskip],cmap=vcmap,pivot="mid",scale=vsize*vscale)
                else:
                    vect = theAxes.quiver(x[::vskip],y[::vskip],u1[::vskip,::vskip],v1[::vskip,::vskip],\
                                          color=vcolor,pivot="mid",scale=vsize*vscale)

                # Plot the scale of the vectors under the axes
                unit_u = self.data[vec+"_"+dir_x]["unit"]
                if vkey:
                    theAxes.quiverkey(vect,vkey_pos[0],vkey_pos[1], vscale,"%.2f [%s]" % (vscale, unit_u),\
                                      labelpos="E", coordinates="axes", color="k", labelcolor="k",zorder=100)

            if stream:
                if scmap:
                    if scmap.startswith("log"):
                        w2 = np.log10(w2)
                        scmap = scmap.split(",")[1]
                    strm = theAxes.streamplot(x,y,u2,v2,color=w2,cmap=scmap)
                else:
                    strm = theAxes.streamplot(x,y,u2,v2,color=scolor)
            
            if self.info["nsinks"] > 0 and sinks:
                sinkMasstot=0.0
                if dz == 0.0:
                    subset = np.where(self.get("r")[cube] < dx*0.01)
                    thickness = 0.5*np.average(celldx[subset])
                else:
                    thickness = 0.5*dz
                div = np.sqrt(a_plane**2 + b_plane**2 + c_plane**2)
                for key in self.sinks.keys():
                    dist = (a_plane*self.sinks[key]["x"]+b_plane*self.sinks[key]["y"]+c_plane*self.sinks[key]["z"]+d_plane) / div
                    if abs(dist) <= thickness:
                        crad = max(self.sinks[key]["radius"],dx*0.01)
                        circle1 = plt.Circle((self.sinks[key][dir_x],self.sinks[key][dir_y]),crad,edgecolor="none",facecolor="w",alpha=0.5)
                        circle2 = plt.Circle((self.sinks[key][dir_x],self.sinks[key][dir_y]),crad,facecolor="none",edgecolor="k")
                        circle3 = plt.Circle((self.sinks[key][dir_x],self.sinks[key][dir_y]),crad*0.2,color="k")
                        theAxes.add_patch(circle1)
                        theAxes.add_patch(circle2)
                        theAxes.add_patch(circle3)
                        sinkMasstot+=self.sinks[key]["mass"]
                theAxes.text(0.02,-0.09,"Msink = %4.1f Msun" % sinkMasstot,transform=theAxes.transAxes,color="k")

            try:
                title += ""
                theAxes.set_title(title)
            except TypeError:
                theAxes.set_title("Time = %.3f %s" % (self.info["time"]/conf.constants[conf.default_values["time_unit"]],conf.default_values["time_unit"]))
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
        
    #=======================================================================================
    # Plot a 1D profile through the data cube.
    #=======================================================================================
    def plot_profile(self,direction,var,x=0.0,y=0.0,z=0.0,var_z=None,fname=None,logz=False,\
                     axes=None,copy=False,xmin=None,xmax=None,ymin=None,ymax=None,\
                     new_window=False,update=None,marker=None,iskip=1,cbar=True,\
                     clear=True,plot=True,**kwargs):
        
        # Possibility of updating the data from inside the plotting routines
        try:
            update += 0
            self.update_values(nout=update)
        except TypeError:
            pass
                
        # Get the data values and units
        datax  = self.get(direction)
        datay  = self.get(var)
        xlabel = self.data[direction]["label"]+" ["+self.data[direction]["unit"]+"]"
        ylabel = self.data[var]["label"]+" ["+self.data[var]["unit"]+"]"
        #if var_z:
            #dataz  = self.data[var_z]["values"]
            #zlabel = self.data[var_z]["label"]
        
        # Define plotting range
        autoxmin = False
        autoxmax = False
        autoymin = False
        autoymax = False
        
        try:
            xmin += 0
        except TypeError:
            xmin = np.amin(datax)
            autoxmin = True
        try:
            xmax += 0
        except TypeError:
            xmax = np.amax(datax)
            autoxmax = True
        try:
            ymin += 0
        except TypeError:
            ymin = np.amin(datay)
            autoymin = True
        try:
            ymax += 0
        except TypeError:
            ymax = np.amax(datay)
            autoymax = True
        
        # Select only the cells in contact with the profile line
        dirs = "xyz".replace(direction,"")
        cube = np.where(np.logical_and(datax-0.5*self.get("dx") <= xmax,\
                        np.logical_and(datax+0.5*self.get("dx") >= xmin,\
                        np.logical_and(abs(self.get(dirs[0])-eval(dirs[0])) <= 0.51*self.get("dx"),\
                                       abs(self.get(dirs[1])-eval(dirs[1])) <= 0.51*self.get("dx")))))
        order = datax[cube].argsort()
        x = datax[cube][order]
        y = datay[cube][order]
        if var_z:
            z  = self.get(var_z)[cube][order]
            zlabel = self.data[var_z]["label"]
        
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
            
            theAxes.plot(x,y,marker=marker,**kwargs)
                            
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
                plt.show(block=False)

        if copy:
            return x,y,z
        else:
            return


