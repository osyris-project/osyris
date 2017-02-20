from pylab import *
import numpy as np

# Read arguments
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = 1

#======================================================================================

# Physical constants
c       = 2.99792458e+10
g       = 6.67259850e-08
hplanck = 6.62607554e-27
#kb      = 1.38065812e-16
mh      = 1.6600000e-24
kb      = 1.380658e-16
a_r     = (8.0*(pi**5)*(kb**4))/(15.0*(hplanck**3)*(c**3))
au      = 1.495980e+13
msol    = 1.989100e+33
lsol    = 3.826800e+33
yr      = 3.155760e+07
mu      = 2.31

#======================================================================================

key = [" "] * 14
key[ 0] = "ilevel"
key[ 1] =  "x"
key[ 2] =  "y"
key[ 3] =  "z"
key[ 4] =  "dx"
key[ 5] =  "rho"
key[ 6] =  "u"
key[ 7] =  "v"
key[ 8] =  "w"
key[ 9] =  "T"
key[10] =  "Bx"
key[11] =  "By"
key[12] =  "Bz"
key[13] =  "B"
nkeys = shape(key)[0]

#======================================================================================

def plot_histogram(data_array,var_x,var_y,fname="fig.pdf",zlog=True):

    no_x = True
    no_y = True
    
    for i in range(nkeys):
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
        
        xmin = amin(data_array[:,ix])
        xmax = amax(data_array[:,ix])
        ymin = amin(data_array[:,iy])
        ymax = amax(data_array[:,iy])

        xe = linspace(xmin,xmax,nx)
        ye = linspace(ymin,ymax,ny)

        z1, yedges1, xedges1 = histogram2d(data_array[:,iy],data_array[:,ix],bins=(ye,xe))

        x = zeros([nx-1])
        y = zeros([ny-1])

        for i in range(nx-1):
            x[i] = 0.5*(xe[i]+xe[i+1])
        for j in range(ny-1):
            y[j] = 0.5*(ye[j]+ye[j+1])

        if zlog:
            z = log10(z1)
            
        fig = matplotlib.pyplot.figure()
        ax  = fig.add_subplot(111)
        ratio = 0.7
        sizex = 10.0
        fig.set_size_inches(sizex,ratio*sizex)
        cont = ax.contourf(x,y,z,20,cmap="YlGnBu")
        ax.set_xlabel(var_x)
        ax.set_ylabel(var_y)
        fig.savefig(fname,bbox_inches="tight")

    return

#======================================================================================

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
        
        #w = 

        z, yedges1, xedges1 = histogram2d(data_array[:,iy],data_array[:,ix],bins=(ye,xe),normed=True,weights=data_array[:,iz])
        
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

#======================================================================================


#======================================================================================


# Load data
sout = str(nout).zfill(5)
datafile = "output_"+sout+"/output_"+sout+".dat"
print("Reading file: "+datafile) 
data = loadtxt(datafile)
#print "converting to log"
data_log = log10(data)
data[:,1] = data[:,1]/au
data[:,2] = data[:,2]/au
data[:,3] = data[:,3]/au

#======================================================================================


# Density - B field
plot_histogram(data_log,"rho","B",fname="brho.pdf")

# Density - Temperature
plot_histogram(data_log,"rho","T",fname="trho.pdf")

# x,y density slice
plot_slice(data,"x","y","rho",fname="rhoxy.pdf",zlog=True,xmin=-50.0,xmax=50.0,ymin=-50.0,ymax=50.0)
