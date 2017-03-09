from pylab import *
import read_ramses_data as rd

#======================================================================================

nvar_lin = 15
nvar_log = 4
nvar_tot = nvar_lin + nvar_log
key = [" "] * nvar_tot
key[ 0] = "ilevel"
key[ 1] =  "x"
key[ 2] =  "y"
key[ 3] =  "z"
key[ 4] =  "dx"
key[ 5] =  "rho"
key[ 6] =  "vel"
key[ 7] =  "T"
key[ 8] =  "B"
key[ 9] =  "u"
key[10] =  "v"
key[11] =  "w"
key[12] =  "Bx"
key[13] =  "By"
key[14] =  "Bz"
key[15] =  "logrho"
key[16] =  "logvel"
key[17] =  "logT"
key[18] =  "logB"

#======================================================================================

def ramses_output(nout=1,xc=0.5,yc=0.5,zc=0.5,lmax=0):
    infile = "output_"+str(nout).zfill(5)
    [data1,nn] = rd.ramses_data(infile,xc,yc,zc,lmax)
    data2 = data1[:nn,:]
    return data2

#======================================================================================

def plot_histogram(data_array,var_x,var_y,fname=None,zlog=True,axes=None):

    no_x = True
    no_y = True
    
    for i in range(nvar_tot):
        if var_x == key[i]:
            no_x = False
            ix = i
        if var_y == key[i]:
            no_y = False
            iy = i
    
    if no_x:
        print("Bad value for x axis")
    elif no_y:
        print("Bad value for y axis")
    else:
        
        # Parameters
        nx = 101
        ny = 101
        
        if ix >= nvar_lin:
            ix = ix-10
            datax = log10(data_array[:,ix])
        else:
            datax = data_array[:,ix]
        if iy >= nvar_lin:
            iy = iy-10
            datay = log10(data_array[:,iy])
        else:
            datay = data_array[:,iy]
        
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
            cont = axes.contourf(x,y,z,20,cmap="YlGnBu")
            axes.set_xlabel(var_x)
            axes.set_ylabel(var_y)
        else:
            fig = matplotlib.pyplot.figure()
            ax  = fig.add_subplot(111)
            ratio = 0.7
            sizex = 10.0
            fig.set_size_inches(sizex,ratio*sizex)
            cont = ax.contourf(x,y,z,20,cmap="YlGnBu")
            ax.set_xlabel(var_x)
            ax.set_ylabel(var_y)
            if fname:
                fig.savefig(fname,bbox_inches="tight")
            else:
                show()

    return

#======================================================================================

def plot_slice(data_array,var_x,var_y,var_z,fname="fig.pdf",zlog=True,xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0):

    no_x = True
    no_y = True
    no_z = True
    
    for i in range(nkeys):
        if var_x == key[i]:
            no_x = False
            ix = i
        if var_y == key[i]:
            no_y = False
            iy = i
        if var_z == key[i]:
            no_z = False
            iz = i
    
    if no_x:
        print("Bad value for x axis: "+var_x)
    elif no_y:
        print("Bad value for y axis: "+var_y)
    elif no_z:
        print("Bad value for z axis: "+var_z)
    else:
        
        # Parameters
        nx = 101
        ny = 101
        
        #xmin = amin(data_array[:,ix])
        #xmax = amax(data_array[:,ix])
        #ymin = amin(data_array[:,iy])
        #ymax = amax(data_array[:,iy])

        xe = linspace(xmin,xmax,nx)
        ye = linspace(ymin,ymax,ny)
        
        data_slice = data_array[where(abs(data_array[:,3]) < 10.0)]

        z, yedges1, xedges1 = histogram2d(data_slice[:,iy],data_slice[:,ix],bins=(ye,xe),normed=True,weights=data_slice[:,iz])
        
        print(amin(z),amax(z))

        x = zeros([nx-1])
        y = zeros([ny-1])

        for i in range(nx-1):
            x[i] = 0.5*(xe[i]+xe[i+1])
        for j in range(ny-1):
            y[j] = 0.5*(ye[j]+ye[j+1])

        if zlog:
            z = log10(z)
            
        fig = matplotlib.pyplot.figure()
        ax  = fig.add_subplot(111)
        #ratio = 0.7
        #sizex = 10.0
        #fig.set_size_inches(sizex,ratio*sizex)
        cont = ax.contourf(x,y,z,20,cmap="YlGnBu")
        ax.set_xlabel(var_x)
        ax.set_ylabel(var_y)
        cbar = fig.colorbar(cont)
        cbar.ax.set_ylabel(var_z)

        fig.savefig(fname,bbox_inches="tight")

    return


    
