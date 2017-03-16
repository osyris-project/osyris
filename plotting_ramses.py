import numpy as np
import glob
import read_ramses_data as rd
import matplotlib.pyplot as plt

#=======================================================================================

class RamsesOutput:
 
    def __init__(self,nout=1,maxlevel=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale="cm"):
        
        if nout == -1:
            filelist = sorted(glob.glob("output*"))
            infile = filelist[-1]
        else:
            infile = "output_"+str(nout).zfill(5)
        
        try:
            center += 0
            xc = center[0]
            yc = center[1]
            zc = center[2]
        except TypeError:
            xc = 0.5
            yc = 0.5
            zc = 0.5
        
        scalelist = {"cm": 1.0, "au": 1.495980e+13, "pc": 3.085678e+18}
        
        [data1,names,nn,ncpu,ndim,lmin,lmax,nstep,boxsize,time,ud,ul,ut] = rd.ramses_data(infile,maxlevel,xc,yc,zc,dx,dy,dz,scalelist[scale])
        
        self.data = dict()
        
        self.data["info"] = dict()
        self.data["info"]["ncells"]  = nn
        self.data["info"]["ncpu"]    = ncpu
        self.data["info"]["ndim"]    = ndim
        self.data["info"]["lmin"]    = lmin
        self.data["info"]["lmax"]    = lmax
        self.data["info"]["nstep"]   = nstep
        self.data["info"]["boxsize"] = boxsize
        self.data["info"]["time"]    = time
        self.data["info"]["ud"]      = ud
        self.data["info"]["ul"]      = ul
        self.data["info"]["ut"]      = ut
        
        list_vars = names.split()
        for i in range(len(list_vars)):
            theKey = list_vars[i]
            self.data[theKey] = dict()
            [norm,uu] = get_units(theKey,ud,ul,ut,scale)
            self.data[theKey]["values"] = data1[:nn,i]*norm
            self.data[theKey]["unit"  ] = uu
            self.data[theKey]["label" ] = theKey

        # Modifications for coordinates and cell sizes
        try:
            lc = len(center)
            if center == "auto":
                maxloc = np.argmax(self.data["density"]["values"])
                xc = self.data["x"]["values"][maxloc]
                yc = self.data["y"]["values"][maxloc]
                zc = self.data["z"]["values"][maxloc]
            elif lc == 3:
                xc = center[0]*self.data["info"]["boxsize"]
                yc = center[1]*self.data["info"]["boxsize"]
                zc = center[2]*self.data["info"]["boxsize"]
            else:
                print "Bad center value"
                return
        except TypeError:
            xc = 0.5*self.data["info"]["boxsize"]
            yc = 0.5*self.data["info"]["boxsize"]
            zc = 0.5*self.data["info"]["boxsize"]
        self.data["x"]["values"] = (self.data["x"]["values"] - xc)/scalelist[scale]
        self.data["y"]["values"] = (self.data["y"]["values"] - yc)/scalelist[scale]
        if ndim > 2:
            self.data["z"]["values"] = (self.data["z"]["values"] - zc)/scalelist[scale]
        
        self.data["dx"]["values"] = self.data["dx"]["values"]/scalelist[scale]
        
        # Magnetic field
        self.data["B_x"] = dict()
        self.data["B_x"]["values"] = 0.5*(self.data["B_left_x"]["values"]+self.data["B_right_x"]["values"])
        self.data["B_x"]["unit"  ] = "G"
        self.data["B_x"]["label" ] = "B_x"
        
        self.data["B_y"] = dict()
        self.data["B_y"]["values"] = 0.5*(self.data["B_left_y"]["values"]+self.data["B_right_y"]["values"])
        self.data["B_y"]["unit"  ] = "G"
        self.data["B_y"]["label" ] = "B_y"
        
        self.data["B_z"] = dict()
        self.data["B_z"]["values"] = 0.5*(self.data["B_left_z"]["values"]+self.data["B_right_z"]["values"])
        self.data["B_z"]["unit"  ] = "G"
        self.data["B_z"]["label" ] = "B_z"
        
        self.data["B"] = dict()
        self.data["B"]["values"] = np.sqrt(self.data["B_x"]["values"]**2+self.data["B_y"]["values"]**2+self.data["B_z"]["values"]**2)
        self.data["B"]["unit"  ] = "G"
        self.data["B"]["label" ] = "B"
    
    #------------------------------------------------------------------------------
    
    def get_values(self,variable,key="values"):
        return self.data[variable][key]
    
    #------------------------------------------------------------------------------
    
    def get(self,variable):
        return self.data[variable]
    
    #------------------------------------------------------------------------------
    
    def new_field(self,name,values,unit,label):
        self.data[name] = dict()
        self.data[name]["values"] = values
        self.data[name]["unit"  ] = unit
        self.data[name]["label" ] = label
        
#======================================================================================

def plot_histogram(var_x,var_y,var_z=None,fname=None,logz=True,axes=None,cmap=None):

    # Parameters
    nx = 101
    ny = 101
    
    try:
        datax = var_x["values"]
        xlabel = var_x["label"]+" ["+var_x["unit"]+"]"
    except IndexError:
        datax = var_x
        xlabel = ""
        
    try:
        datay = var_y["values"]
        ylabel = var_y["label"]+" ["+var_y["unit"]+"]"
    except IndexError:
        datay = var_y
        ylabel = ""
            
    xmin = np.amin(datax)
    xmax = np.amax(datax)
    ymin = np.amin(datay)
    ymax = np.amax(datay)

    xe = np.linspace(xmin,xmax,nx)
    ye = np.linspace(ymin,ymax,ny)

    z, yedges1, xedges1 = np.histogram2d(datay,datax,bins=(ye,xe))

    contourz = True
    try:
        dataz = var_z["values"]
        zlabel = True
    except IndexError:
        dataz = var_z
        zlabel = False
    except TypeError:
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
        if zlabel:
            leg = [over.collections[0]]
            theplot.legend(leg,[var_z["label"]],loc=2)
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

#======================================================================================

def plot_slice(ramsesdata,var="rho",direction="z",vec=None,streamlines=False,fname=None,dx=1.0,dy=1.0,cmap=None,axes=None,resolution=128):
    
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

    # Make a guess for slice thickness
    dz = 0.05*(0.5*(dx+dy))
    
    data1 = ramsesdata.get_values(dir_x)
    data2 = ramsesdata.get_values(dir_y)
    data3 = ramsesdata.get_values(direction)
    cube  = np.where(np.logical_and(abs(data1) < 0.5*dx,np.logical_and(abs(data2) < 0.5*dy,abs(data3) < 0.5*dz)))
    datax = data1[cube]
    datay = data2[cube]
    dataz = ramsesdata.get_values(var)[cube]
    if vec:
        datau = ramsesdata.get_values(vec+"_"+dir_x)[cube]
        datav = ramsesdata.get_values(vec+"_"+dir_y)[cube]
    celldx = ramsesdata.get_values("dx")[cube]
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
    
    xlab = ramsesdata.get_values(dir_x,key="label")+" ["+ramsesdata.get_values(dir_x,key="unit")+"]"
    ylab = ramsesdata.get_values(dir_y,key="label")+" ["+ramsesdata.get_values(dir_y,key="unit")+"]"
    zlab = ramsesdata.get_values(var,key="label")+" ["+ramsesdata.get_values(var,key="unit")+"]"
    
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

#======================================================================================

def get_units(string,ud,ul,ut,scale="cm"):
    if string == "density":
        return [ud,"g/cm3"]
    elif string[0:8] == "velocity":
        return [ul/ut,"cm/s"]
    elif string[0:2] == "B_":
        return [np.sqrt(4.0*np.pi*ud*(ul/ut)**2),"G"]
    elif string[0:8] == "thermal_pressure":
        return [ud*((ul/ut)**2),"g/cm/s2"]
    elif string == "x":
        return [ul,scale]
    elif string == "y":
        return [ul,scale]
    elif string == "z":
        return [ul,scale]
    elif string == "dx":
        return [ul,scale]
    else:
        return [1.0,""]
