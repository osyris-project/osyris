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
    elif key=="u":
        return  9
    elif key=="v":
        return 10
    elif key=="w":
        return 11
    elif key=="Bx":
        return 12
    elif key=="By":
        return 13
    elif key=="Bz":
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
    elif key=="u":
        return data_array[:, 9]
    elif key=="v":
        return data_array[:,10]
    elif key=="w":
        return data_array[:,11]
    elif key=="Bx":
        return data_array[:,12]
    elif key=="By":
        return data_array[:,13]
    elif key=="Bz":
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

def plot_slice(data_array,dir_z="z",var="rho",fname=None,dx=1.0,dy=1.0,cmap=None,axes=None):

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
    
    if len(dir_z) > 1:
        dir_x = dir_x+dir_z[1:]
        dir_y = dir_y+dir_z[1:]
    
    data1 = get_variable(data_array,dir_x)
    data2 = get_variable(data_array,dir_y)
    data3 = get_variable(data_array,dir_z)
    
    # Make a guess for slice thickness
    dz = 0.05*(0.5*(dx+dy))
    
    #cube = where(logical_and(abs(data1) < 0.5*dx,abs(data2) < 0.5*dy,abs(data3) < 0.5*dz))
    cube = where(logical_and(abs(data1) < 0.5*dx,logical_and(abs(data2) < 0.5*dy,abs(data3) < 0.5*dz)))
    #cube = where(abs(data3) < 0.5*dz)
    
    datax = get_variable(data_array,dir_x)[cube]
    datay = get_variable(data_array,dir_y)[cube]
    dataz = get_variable(data_array,var  )[cube]
    
    print shape(datax)
    
    xmin = -0.5*dx
    xmax =  0.5*dx
    ymin = -0.5*dy
    ymax =  0.5*dy
    
    # Make a guess for lmin
    lmin = 13
    # Find locations where l < lmin
    datal = get_variable(data_array,'level')[cube]
    #print shape(datal)
    locs = where(logical_or(datal < lmin,datal>=lmin-3))
    print locs,shape(locs)
    [dum,nl] = shape(locs)
    nnew = 0
    for l in range(nl):
        #print 'l is',l,nl
        #print locs[0][1]
        k = locs[0][l]
        #print datal[k]
        nextra = 2**(lmin-int(datal[k]))
        xc = datax[k]
        yc = datay[k]
        vc = dataz[k]
        nc = shape(datax)[0]
        #print shape(datax)[0],nc
        datax = resize(datax,(nc+nextra**2))
        datay = resize(datay,(nc+nextra**2))
        dataz = resize(dataz,(nc+nextra**2))
        #print shape(datax)
        dxcell = data_array[k,get_index("dx")]/au
        dxnew = dxcell/nextra
        ic = 0
        for j in range(nextra):
            yy = yc - 0.5*dxcell + (j+0.5)*dxnew
            for i in range(nextra):
                xx = xc - 0.5*dxcell + (i+0.5)*dxnew
                ic = ic + 1
                datax[nc+ic-1] = xx
                datay[nc+ic-1] = yy
                dataz[nc+ic-1] = vc
                nnew = nnew + 1
    
    print shape(datax)
    print datax[-1],datay[-1],dataz[-1]
    
    #print "added ",nnew," new cells"
    
    ## Parameters
    ## First pass: low resolution background
    #nx = 51
    #ny = 51
    #xe = linspace(xmin,xmax,nx)
    #ye = linspace(ymin,ymax,ny)
    #denominator1, xedges1, yedges1 = histogram2d(datay,datax,bins=[ye, xe])
    #nominator1  , xedges1, yedges1 = histogram2d(datay,datax,bins=[ye, xe], weights=dataz)
    #z1 = nominator1 / denominator1
    #x1 = zeros([nx-1])
    #y1 = zeros([ny-1])
    #for i in range(nx-1):
        #x1[i] = 0.5*(xe[i]+xe[i+1])
    #for j in range(ny-1):
        #y1[j] = 0.5*(ye[j]+ye[j+1])
    # Second pass: high resolution foreground
    nx = 101
    ny = 101
    xe = linspace(xmin,xmax,nx)
    ye = linspace(ymin,ymax,ny)
    denominator, xedges1, yedges1 = histogram2d(datay,datax,bins=[ye, xe])
    nominator  , xedges1, yedges1 = histogram2d(datay,datax,bins=[ye, xe], weights=dataz)
    z = nominator / denominator
    #zmin = amin(z2[where(denominator > 0.0)])
    #zmax = amax(z2[where(denominator > 0.0)])
    [nny,nnx] = shape(denominator)
    x = zeros([nnx])
    y = zeros([nny])
    for i in range(nnx):
        x[i] = 0.5*(xe[i]+xe[i+1])
    for j in range(nny):
        y[j] = 0.5*(ye[j]+ye[j+1])

    ## Fill empty pixels
    #if amin(denominator) < 1.0:
        #print shape(denominator)
        #[nny,nnx] = shape(denominator)
        ## Make several iterations for large empty areas
        #empty_cells_exist = True
        #niter = 0
        #while empty_cells_exist:
            #niter = niter + 1
            #empty_cells_exist = False
            #nfix=0
            #for j in range(nny):
                #for i in range(nnx):
                    #if denominator[j,i] < 1.0:
                        ##print x[i],y[j]
                        #d_loc = 0.0
                        #n_loc = 0.0
                        #count = 0.0
                        ##print max(0,j-1),min(nny-1,j+1),max(0,i-1),min(nnx-1,i+1)
                        #for jj in (max(0,j-1),min(nny-1,j+1)):
                            #for ii in (max(0,i-1),min(nnx-1,i+1)):
                                #if denominator[jj,ii] > 0.0:
                                    #d_loc = d_loc + denominator[jj,ii]
                                    #n_loc = n_loc +   nominator[jj,ii]
                                    #count = count + 1.0
                                    ##print i,j,ii,jj,d_loc,n_loc,count,denominator[jj,ii],nominator[jj,ii]
                        #if count > 0.0:
                            #denominator[j,i] = d_loc / count
                            #nominator[j,i] = n_loc / count
                            #empty_cells_exist = True
                            #nfix = nfix + 1
                        
                        ##input("Press Enter to continue...")
            ##print niter,empty_cells_exist,nfix
            ##break
            
    #z = nominator / denominator





    #print zmin,zmax
    
    #levs = linspace(zmin,zmax,20)
    
    if axes:
        cont = axes.contourf(x,y,z,20,cmap=cmap)
        #cont = axes.imshow(z1,cmap=cmap,interpolation='None',extent=[xmin,xmax,ymin,ymax])
        #cont = axes.imshow(z,cmap=cmap,interpolation='None',extent=[xmin,xmax,ymin,ymax])
        #axes.scatter(datax, datay, c=dataz)
        #cont2 = axes.contourf(x2,y2,z2,levels=levs,cmap=cmap)
        axes.set_xlabel(dir_x)
        axes.set_ylabel(dir_y)
        cbar = colorbar(cont)
    else:
        fig = matplotlib.pyplot.figure()
        ax  = fig.add_subplot(111)
        #ratio = 0.7
        #sizex = 10.0
        #fig.set_size_inches(sizex,ratio*sizex)
        cont = ax.contourf(x,y,z,20,cmap=cmap)
        #cont2 = ax.contourf(x2,y2,z2,levels=levs,cmap=cmap)
        cbar = fig.colorbar(cont)
        ax.set_xlabel(dir_x)
        ax.set_ylabel(dir_y)
        if fname:
            fig.savefig(fname,bbox_inches="tight")
        else:
            show()
        
    return


    
