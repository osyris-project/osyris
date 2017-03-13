from pylab import *
import glob
import read_ramses_data as rd

#=======================================================================================

class RamsesOutput:
 
    def __init__(self,nout=1,maxlevel=0,center=None,dx=0.0,dy=0.0,dz=0.0,scale=1.0):
        
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
            
        [data1,nn,ncpu,ndim,lmin,lmax,nstep,boxsize,time] = rd.ramses_data(infile,maxlevel,xc,yc,zc,dx,dy,dz,scale)
        
        self.ncells  = nn
        self.ncpu    = ncpu
        self.ndim    = ndim
        self.lmin    = lmin
        self.lmax    = lmax
        self.nstep   = nstep
        self.boxsize = boxsize
        self.time    = time
        
        self.level  = data1[:nn, 0]
        self.dx     = data1[:nn, 4]/scale
        self.rho    = data1[:nn, 5]
        self.T      = data1[:nn, 7]
        
        self.x = zeros([nn,4])
        self.x[:,0] = data1[:nn, 1]
        self.x[:,1] = data1[:nn, 2]
        self.x[:,2] = data1[:nn, 3]
        self.x[:,3] = data1[:nn, 4]/scale
        
        if center == "auto":
            maxloc = argmax(self.rho)
            xc = self.x[maxloc,0]
            yc = self.x[maxloc,1]
            zc = self.x[maxloc,2]
        elif len(center) == 3:
            xc = center[0]*self.boxsize
            yc = center[1]*self.boxsize
            zc = center[2]*self.boxsize
        else:
            xc = 0.5*self.boxsize
            yc = 0.5*self.boxsize
            zc = 0.5*self.boxsize
        self.x[:,0] = (self.x[:,0] - xc)/scale
        self.x[:,1] = (self.x[:,1] - yc)/scale
        self.x[:,2] = (self.x[:,2] - zc)/scale
        
        self.vel = zeros([nn,4])
        self.vel[:,0] = data1[:nn, 9]
        self.vel[:,1] = data1[:nn,10]
        self.vel[:,2] = data1[:nn,11]
        self.vel[:,3] = data1[:nn, 6]
        
        self.B = zeros([nn,4])
        self.B[:,0] = data1[:nn,12]
        self.B[:,1] = data1[:nn,13]
        self.B[:,2] = data1[:nn,14]
        self.B[:,3] = data1[:nn, 8]
        
#======================================================================================

def plot_histogram(datax,datay,dataz=None,fname=None,zlog=True,axes=None,cmap=None):

    # Parameters
    nx = 101
    ny = 101
            
    xmin = amin(datax)
    xmax = amax(datax)
    ymin = amin(datay)
    ymax = amax(datay)

    xe = linspace(xmin,xmax,nx)
    ye = linspace(ymin,ymax,ny)

    z, yedges1, xedges1 = histogram2d(datay,datax,bins=(ye,xe))

    try:
        contourz = (len(dataz) > 0)
    except TypeError:
        contourz = False

    if contourz:
        z1, yedges1, xedges1 = histogram2d(datay,datax,bins=(ye,xe),weights=dataz)
        z2 = z1/z
        zmin = amin(dataz)
        zmax = amax(dataz)

    x = zeros([nx-1])
    y = zeros([ny-1])

    for i in range(nx-1):
        x[i] = 0.5*(xe[i]+xe[i+1])
    for j in range(ny-1):
        y[j] = 0.5*(ye[j]+ye[j+1])

    if zlog:
        z = log10(z)
    
    if axes:
        ax = axes
    else:
        fig = matplotlib.pyplot.figure()
        ax  = fig.add_subplot(111)
    
    cont = ax.contourf(x,y,z,20,cmap=cmap)
    if contourz:
        over = ax.contour(x,y,z2,levels=arange(zmin,zmax+1),colors='k')
        ax.clabel(over,inline=1,fmt='%i')
    
    if fname:
        fig.savefig(fname,bbox_inches="tight")
    elif axes:
        pass
    else:
        show()

    return

#======================================================================================

def plot_slice(xyz,var,direction=2,vec=None,streamlines=False,fname=None,dx=1.0,dy=1.0,cmap=None,axes=None,resolution=128):
    
    try:
        vectors = (len(vec) > 0)
    except TypeError:
        vectors = False

    # Make a guess for slice thickness
    dz = 0.05*(0.5*(dx+dy))
    
    dirs = [(direction+1)%3,(direction+2)%3]
    ix = amin(dirs)
    iy = amax(dirs)
    
    data1 = xyz[:,ix]
    data2 = xyz[:,iy]
    data3 = xyz[:,direction]
    cube  = where(logical_and(abs(data1) < 0.5*dx,logical_and(abs(data2) < 0.5*dy,abs(data3) < 0.5*dz)))
    datax = data1[cube]
    datay = data2[cube]
    dataz = var  [cube]
    if vectors:
        datau  = vec[:,ix][cube]
        datav  = vec[:,iy][cube]
    celldx = xyz[:,3][cube]
    ncells = shape(datax)[0]
    
    xmin = -0.5*dx
    xmax =  0.5*dx
    ymin = -0.5*dy
    ymax =  0.5*dy
    
    nx = resolution
    ny = resolution
    dpx = (xmax-xmin)/nx
    dpy = (ymax-ymin)/ny
    
    z1 = zeros([ny,nx])
    z2 = zeros([ny,nx])
    if vectors:
        u1 = zeros([ny,nx])
        v1 = zeros([ny,nx])
        z3 = zeros([ny,nx])
    
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
                if vectors:
                    u1[j,i] = u1[j,i] + datau[n]
                    v1[j,i] = v1[j,i] + datav[n]
                    z3[j,i] = z3[j,i] + sqrt(datau[n]**2+datav[n]**2)
        
    z = z1/z2
    if vectors:
        u = u1/z2
        v = v1/z2
        w = z3/z2
    
    x = linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
    y = linspace(ymin+0.5*dpy,ymax-0.5*dpy,ny)
    iskip = int(0.071*resolution)
    
    if axes:
        ax = axes
    else:
        fig = matplotlib.pyplot.figure()
        ax  = fig.add_subplot(111)
        
    cont = ax.contourf(x,y,z,20,cmap=cmap)
    cbar = colorbar(cont,ax=ax)
    if vectors:
        if streamlines:
            if streamlines == "log":
                w = log10(w)
            strm = ax.streamplot(x,y,u,v,color=w,cmap='Greys')
        else:
            vect = ax.quiver(x[::iskip],y[::iskip],u[::iskip,::iskip],v[::iskip,::iskip],w[::iskip,::iskip],cmap='Greys',pivot='mid')
    #ax.set_xlabel(dir_x)
    #ax.set_ylabel(dir_y)
    if fname:
        fig.savefig(fname,bbox_inches="tight")
    elif axes:
        pass
    else:
        show()
        
    return
