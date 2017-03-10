from pylab import *
import glob
import read_ramses_data as rd

# === Physical constants ==============================================================

au = 1.495980e+13

##======================================================================================

def get_index(key="rho"):
    if key=="level":
        return 0
    elif key=="x":
        return  1
    elif key=="y":
        return  2
    elif key=="z":
        return  3
    elif key=="dx":
        return  4
    elif key=="rho":
        return  5
    elif key=="vel":
        return  6
    elif key=="T":
        return  7
    elif key=="B":
        return  8
    elif key=="vel_x":
        return  9
    elif key=="vel_y":
        return 10
    elif key=="vel_z":
        return 11
    elif key=="B_x":
        return 12
    elif key=="B_y":
        return 13
    elif key=="B_z":
        return 14
    else:
       print("Unknown key")

#======================================================================================

def get_variable(data_array,key="rho"):
    if key=="level":
        return data_array[:, 0]
    elif key=="x":
        return data_array[:, 1]
    elif key=="y":
        return data_array[:, 2]
    elif key=="z":
        return data_array[:, 3]
    elif key=="x_au":
        return data_array[:, 1]/au
    elif key=="y_au":
        return data_array[:, 2]/au
    elif key=="z_au":
        return data_array[:, 3]/au
    elif key=="dx":
        return data_array[:, 4]
    elif key=="dx_au":
        return data_array[:, 4]/au
    elif key=="rho":
        return data_array[:, 5]
    elif key=="vel":
        return data_array[:, 6]
    elif key=="T":
        return data_array[:, 7]
    elif key=="B":
        return data_array[:, 8]
    elif key=="vel_x":
        return data_array[:, 9]
    elif key=="vel_y":
        return data_array[:,10]
    elif key=="vel_z":
        return data_array[:,11]
    elif key=="B_x":
        return data_array[:,12]
    elif key=="B_y":
        return data_array[:,13]
    elif key=="B_z":
        return data_array[:,14]
    elif key=="logrho":
        return log10(data_array[:,5])
    elif key=="logvel":
        return log10(data_array[:,6])
    elif key=="logT":
        return log10(data_array[:,7])
    elif key=="logB":
        return log10(data_array[:,8])
    else:
        print("Cannot fetch variable, unrecognised data key: "+key)

#======================================================================================

def ramses_output(nout=1,xc=0.5,yc=0.5,zc=0.5,lmax=0):
    if nout == -1:
        filelist = sorted(glob.glob("output*"))
        infile = filelist[-1]
    else:
        infile = "output_"+str(nout).zfill(5)
    [data1,nn] = rd.ramses_data(infile,xc,yc,zc,lmax)
    data2 = data1[:nn,:]
    return data2

#======================================================================================

def plot_histogram(data_array,var_x,var_y,fname=None,zlog=True,axes=None,cmap=None):

    # Parameters
    nx = 101
    ny = 101
    
    datax = get_variable(data_array,var_x)
    datay = get_variable(data_array,var_y)
        
    xmin = amin(datax)
    xmax = amax(datax)
    ymin = amin(datay)
    ymax = amax(datay)

    xe = linspace(xmin,xmax,nx)
    ye = linspace(ymin,ymax,ny)

    z, yedges1, xedges1 = histogram2d(datay,datax,bins=(ye,xe))

    x = zeros([nx-1])
    y = zeros([ny-1])

    for i in range(nx-1):
        x[i] = 0.5*(xe[i]+xe[i+1])
    for j in range(ny-1):
        y[j] = 0.5*(ye[j]+ye[j+1])

    if zlog:
        z = log10(z)
    
    if axes:
        cont = axes.contourf(x,y,z,20,cmap=cmap)
        axes.set_xlabel(var_x)
        axes.set_ylabel(var_y)
    else:
        fig = matplotlib.pyplot.figure()
        ax  = fig.add_subplot(111)
        ratio = 0.7
        sizex = 10.0
        fig.set_size_inches(sizex,ratio*sizex)
        cont = ax.contourf(x,y,z,20,cmap=cmap)
        ax.set_xlabel(var_x)
        ax.set_ylabel(var_y)
        if fname:
            fig.savefig(fname,bbox_inches="tight")
        else:
            show()

    return

#======================================================================================

def plot_slice(data_array,dir_z="z",var="rho",vec=None,streamlines=False,fname=None,dx=1.0,dy=1.0,cmap=None,axes=None,resolution=128):

    if dir_z[0]=="z":
        dir_x = "x"
        dir_y = "y"
    elif dir_z[0]=="x":
        dir_x = "y"
        dir_y = "z"
    elif dir_z[0]=="y":
        dir_x = "x"
        dir_y = "z"
    else:
        print("Bad z direction: "+dir_z)
    
    if vec:
        vec_x = vec+"_"+dir_x
        vec_y = vec+"_"+dir_y
        
    #if stream:
        #str_x = stream+"_"+dir_x
        #str_y = stream+"_"+dir_y
    
    
    if len(dir_z) > 1:
        dir_x = dir_x+dir_z[1:]
        dir_y = dir_y+dir_z[1:]
    
    data1 = get_variable(data_array,dir_x)
    data2 = get_variable(data_array,dir_y)
    data3 = get_variable(data_array,dir_z)
    
    # Make a guess for slice thickness
    dz = 0.05*(0.5*(dx+dy))
    
    cube = where(logical_and(abs(data1) < 0.5*dx,logical_and(abs(data2) < 0.5*dy,abs(data3) < 0.5*dz)))
    
    datax  = get_variable(data_array,dir_x)[cube]
    datay  = get_variable(data_array,dir_y)[cube]
    dataz  = get_variable(data_array,var  )[cube]
    if vec:
        datau  = get_variable(data_array,vec_x)[cube]
        datav  = get_variable(data_array,vec_y)[cube]
    celldx = get_variable(data_array,"dx_au")[cube]
    
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
    if vec:
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
                if vec:
                    u1[j,i] = u1[j,i] + datau[n]
                    v1[j,i] = v1[j,i] + datav[n]
                    z3[j,i] = z3[j,i] + sqrt(datau[n]**2+datav[n]**2)
                
    z = z1/z2
    if vec:
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
    if vec:
        if streamlines:
            if streamlines == "log":
                w = log10(w)
            strm = ax.streamplot(x,y,u,v,color=w,cmap='Greys')
        else:
            vect = ax.quiver(x[::iskip],y[::iskip],u[::iskip,::iskip],v[::iskip,::iskip],w[::iskip,::iskip],cmap='Greys',pivot='mid')
    ax.set_xlabel(dir_x)
    ax.set_ylabel(dir_y)
    if fname:
        fig.savefig(fname,bbox_inches="tight")
    elif axes:
        pass
    else:
        show()
        
    return


    
