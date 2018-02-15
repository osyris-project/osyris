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

from load_ramses_data import get_binary_data
import config_osiris as conf
import numpy as np
import struct
from scipy.interpolate import RegularGridInterpolator

#===================================================================================
# An empty container class which will be filled up as the EOS, opacity, resistivity
# tables are read.
#===================================================================================
class IsmTable():
    
    def __init__(self):
                
        return


#===================================================================================
# Generic function to interpolate simulation data onto opacity table
#===================================================================================
def ism_interpolate(table_container=None,values=[0],points=[0],in_log=False):
    
    func = RegularGridInterpolator(table_container.grid,values)
    
    if in_log:
        return func(points)
    else:
        return np.power(10.0,func(points))
    


#===================================================================================
# Function to read in binary file containing EOS table
#===================================================================================
def read_eos_table(fname="tab_eos.dat"):
    
    print("Loading EOS table: "+fname)
    
    # Read binary EOS file
    with open(fname, mode='rb') as eos_file:
        eosContent = eos_file.read()
    eos_file.close()
    
    # Define data fields. Note that the order is important!
    data_fields = ["rho_eos","ener_eos","temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]

    # Create table container
    theTable = IsmTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    # Get table dimensions
    [theTable.nRho,theTable.nEner] = get_binary_data(fmt="2i",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # Get table limits
    ninteg += 2
    nlines += 1
    [theTable.rhomin,theTable.rhomax,theTable.emin,theTable.emax,theTable.yHe] = \
        get_binary_data(fmt="5d",content=eosContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
        
    array_size = theTable.nRho*theTable.nEner
    array_fmt  = "%id" % array_size
    nfloat += 5
    nlines += 1
    
    # Now loop through all the data fields
    for i in range(len(data_fields)):
        setattr(theTable,data_fields[i],np.reshape(get_binary_data(fmt=array_fmt,content=eosContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),[theTable.nRho,theTable.nEner],order="F"))
        nfloat += array_size
        nlines += 1
    
    del eosContent
    
    theTable.grid = (np.log10(theTable.rho_eos[:,0]), np.log10(theTable.ener_eos[0,:]/theTable.rho_eos[0,:]))
        
    print("EOS table read successfully")

    return theTable


#===================================================================================
# Function to interpolate simulation data onto EOS table
#===================================================================================
def get_eos(holder,eos_fname="tab_eos.dat",variables=["temp_eos","pres_eos","s_eos","cs_eos","xH_eos","xH2_eos","xHe_eos","xHep_eos"]):
    
    if holder.info["eos"] == 0:
        print("Simulation data did not use a tabulated EOS. Exiting")
        return
    else:
        try:
            n = holder.eos_table.nRho
        except AttributeError:
            holder.eos_table = read_eos_table(fname=eos_fname)

        pts = np.array([np.log10(holder.density.values), np.log10(holder.internal_energy.values)]).T
        for var in variables:
            print("Interpolating "+var)
            vals = ism_interpolate(holder.eos_table,np.log10(getattr(holder.eos_table,var)),pts)
            holder.new_field(name=var,label=var,values=vals,verbose=False)

    return







#===================================================================================
#===================================================================================
#===================================================================================

#===================================================================================
# Function to read in binary file containing opacity table
#===================================================================================
def read_opacity_table(fname="vaytet_grey_opacities3D.bin"):
    
    print("Loading opacity table: "+fname)
    
    # Read binary resistivity file
    with open(fname, mode='rb') as kappa_file:
        kappaContent = kappa_file.read()
    kappa_file.close()
    
    # Create table container
    theTable = IsmTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    # Get table dimensions
    [theTable.nx,theTable.ny,theTable.nz] = get_binary_data(fmt="3i",content=kappaContent)
    
    # Get table limits
    [theTable.dx,theTable.dy,theTable.dz,theTable.xmin,theTable.xmax,theTable.ymin,theTable.ymax,theTable.zmin,theTable.zmax] = \
        get_binary_data(fmt="9d",content=kappaContent,correction=12)
    
    # Read table coordinates:
    
    # x: density
    ninteg += 3
    nfloat += 9
    nlines += 1
    theTable.dens = get_binary_data(fmt="%id"%theTable.nx,content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # y: gas temperature
    nfloat += theTable.nx
    nlines += 1
    theTable.tgas = get_binary_data(fmt="%id"%theTable.ny,content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # z: radiation temperature
    nfloat += theTable.ny
    nlines += 1
    theTable.trad = get_binary_data(fmt="%id"%theTable.nz,content=kappaContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # Now read opacities
    array_size = theTable.nx*theTable.ny*theTable.nz
    array_fmt  = "%id" % array_size
    
    #print theTable.nx,theTable.ny,theTable.nz
    
    # Planck mean
    nfloat += theTable.nz
    nlines += 1
    theTable.kappa_p = np.reshape(get_binary_data(fmt=array_fmt,content=kappaContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),[theTable.nx,theTable.ny,theTable.nz],order="F")
    
    # Rosseland mean
    nfloat += array_size
    nlines += 1
    theTable.kappa_r = np.reshape(get_binary_data(fmt=array_fmt,content=kappaContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),[theTable.nx,theTable.ny,theTable.nz],order="F")
    
    del kappaContent
    
    theTable.grid = (theTable.dens,theTable.tgas,theTable.trad)
    
    #print theTable.dens
    #print np.shape(theTable.dens)
    ##print theTable.tgas
    ##print theTable.trad
    
    
    print("Opacity table read successfully")
    
    return theTable


#===================================================================================
# Function to interpolate simulation data onto opacity table
#===================================================================================
def get_opacities(holder,opacity_fname="vaytet_grey_opacities3D.bin",variables=["kappa_p","kappa_r"]):
    
    try:
        n = holder.opacity_table.nx
    except AttributeError:
        holder.opacity_table = read_opacity_table(fname=opacity_fname)
    
    if not hasattr(holder,"radiative_temperature"):
        print("Radiative temperature is not defined. Computing it now.")
        holder.new_field(name="radiative_temperature",operation="(radiative_energy_1/"+str(conf.constants["a_r"])+")**0.25",label="Trad")

    pts = np.array([np.log10(holder.density.values),np.log10(holder.temperature.values),np.log10(holder.radiative_temperature.values)]).T
    for var in variables:
        print("Interpolating "+var)
        vals = ism_interpolate(holder.opacity_table,getattr(holder.opacity_table,var),pts)
        holder.new_field(name=var,label=var,values=vals,verbose=False)

    return



#===================================================================================
#===================================================================================
#===================================================================================

#===================================================================================
# Function to read in binary file containing resistivity table
#===================================================================================
def read_resistivity_table(fname="resistivities_masson2016.bin"):
    
    print("Loading resistivity table: "+fname)
    
    # Read binary resistivity file
    with open(fname, mode='rb') as res_file:
        resContent = res_file.read()
    res_file.close()
    
    # Create table container
    theTable = IsmTable()
    
    # Initialise offset counters and start reading data
    ninteg = nfloat = nlines = nstrin = nquadr = 0
    
    # Get length of record on first line to determine number of dimensions in table
    rec_size = get_binary_data(fmt="i",content=resContent,correction=-4)
    theTable.ndims = rec_size[0]/4
    
    # Get table dimensions
    theTable.nx = np.array(get_binary_data(fmt="%ii"%theTable.ndims,content=resContent))
    #print theTable.nx
    
    # Read table coordinates:
    
    # 1: density
    ninteg += theTable.ndims
    nlines += 1
    theTable.dens = get_binary_data(fmt="%id"%theTable.nx[0],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # 2: gas temperature
    nfloat += theTable.nx[0]
    nlines += 1
    theTable.tgas = get_binary_data(fmt="%id"%theTable.nx[1],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    if theTable.ndims == 4:
        # 3: ionisation rate
        nfloat += theTable.nx[1]
        nlines += 1
        theTable.ionx = get_binary_data(fmt="%id"%theTable.nx[2],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # 4: magnetic field
    nfloat += theTable.nx[-2]
    nlines += 1
    theTable.bmag = get_binary_data(fmt="%id"%theTable.nx[-1],content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat)
    
    # Now read resistivities
    array_size = np.prod(theTable.nx)
    array_fmt  = "%id" % array_size
    
    # Ohmic resistivity
    nfloat += theTable.nx[-1]
    nlines += 1
    theTable.eta_ohm = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")
    
    # Ambipolar resistivity
    nfloat += array_size
    nlines += 1
    theTable.eta_ad = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")
    
    # Hall resistivity
    nfloat += array_size
    nlines += 1
    theTable.eta_hall = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")
    
    # Hall sign
    nfloat += array_size
    nlines += 1
    theTable.eta_hsig = np.reshape(get_binary_data(fmt=array_fmt,content=resContent, \
                ninteg=ninteg,nlines=nlines,nfloat=nfloat),theTable.nx,order="F")
    
    del resContent
    
    if theTable.ndims == 4:
        theTable.grid = (theTable.dens,theTable.tgas,theTable.ionx,theTable.bmag)
    elif theTable.ndims == 3:
        theTable.grid = (theTable.dens,theTable.tgas,theTable.bmag)
    
    # Additional parameters
    theTable.scale_dens = 0.844*2.0/1.66e-24 # 2.0*H2_fraction/mH
    theTable.ionis_rate = 1.0e-17
    
    print("Resistivity table read successfully")
    
    return theTable


#===================================================================================
# Function to interpolate simulation data onto resistivity table
#===================================================================================
def get_resistivities(holder,resistivity_fname="resistivities_masson2016.bin",variables=["eta_ohm","eta_ad","eta_hall"]):
    
    try:
        n = holder.resistivity_table.nx
    except AttributeError:
        holder.resistivity_table = read_resistivity_table(fname=resistivity_fname)
    
    try:
        rho_to_nH = holder.resistivity_table.scale_dens/holder.info["mu_gas"]
    except KeyError:
        rho_to_nH = holder.resistivity_table.scale_dens/2.31 # Default mean atomic weight
    
    if holder.resistivity_table.ndims == 4:
        if not hasattr(holder,"ionisation_rate"):
            holder.new_field(name="ionisation_rate",values=[holder.resistivity_table.ionis_rate]*holder.info["ncells"],label="Ionisation rate")
        pts = np.array([np.log10(holder.density.values*rho_to_nH),np.log10(holder.temperature.values), \
                        np.log10(holder.ionisation_rate.values),np.log10(holder.B.values)]).T
    elif holder.resistivity_table.ndims == 3:
        pts = np.array([np.log10(holder.density.values*rho_to_nH),np.log10(holder.temperature.values), \
                        np.log10(holder.B.values)]).T
    
    for var in variables:
        print("Interpolating "+var)
        vals = ism_interpolate(holder.resistivity_table,getattr(holder.resistivity_table,var),pts)
        if var == "eta_hall":
            hall_sign = np.sign(ism_interpolate(holder.resistivity_table,holder.resistivity_table.eta_hsig,pts,in_log=True))
            holder.new_field(name=var,label=var,values=vals*hall_sign,verbose=False)
        else:
            holder.new_field(name=var,label=var,values=vals,verbose=False)

    return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #print "ndims",ndims
    
    ## Get table dimensions
    #theTable.nx = np.roll(np.array(get_binary_data(fmt="%ii"%ndims,content=resContent)),1)
    #print theTable.nx
    
    #nx_read = np.copy(theTable.nx)
    
    #if ndims == 3:
        #nx_read[0] += 2
    #elif ndims == 4:
        #nx_read[0] += 7
    
    ## Now read the bulk of the table containing abundances in one go
    #ninteg += ndims
    #nlines += 1
    #resistivite_chimie_x = np.reshape(get_binary_data(fmt="%id"%(np.prod(nx_read)),content=resContent,ninteg=ninteg,nlines=nlines,nfloat=nfloat),nx_read,order="F")
    
    #print "kk",theTable.nx,nx_read
    #print resistivite_chimie_x[:,0,0]
    #print np.shape(resistivite_chimie_x)
    
    ## Now we compute conductivities ==========================
    
    ## Define some constants ==================================
    
    ##pi=3.1415927410125732422d0
    #rg      = 0.1e-4       # grain radius 
    #mp      = 1.6726e-24   # proton mass
    #me      = 9.1094e-28   # electron mass
    #mg      = 1.2566e-14   # grain mass
    #e       = 4.803204e-10 # electron charge
    #mol_ion = 29.0*mp      # molecular ion mass
    #Met_ion = 23.5*mp      # atomic ion mass
    #kb      = 1.3807e-16   # Boltzmann
    #clight  = 2.9979250e+10
    ##real(kind=8), parameter ::mH      = 1.6600000d-24
    ##real(kind=8), parameter ::mu_gas = 2.31d0
    ##real(kind=8) :: scale_d = mu_gas*mH
    #rho_s      = 2.3
    #rho_n_tot  = 1.17e-21
    #a_0        = 0.0375e-4
    #a_min      = 0.0181e-4
    #a_max      = 0.9049e-4
    #zeta       = a_min/a_max
    #lambda_pow = -3.5
    
    ## Compute grain distribution =============================
    
    #if ndims == 3:
        #nbins_grains = (theTable.nx[0]-7)/3
        #nion = 7
    #elif ndims == 4:
        #nbins_grains = (theTable.nx[0]-9)/3
        #nion = 9
    #print 'nbins_grains',nbins_grains
    
    #r_g = np.zeros([nbins_grains])
       
    #Lp1 = lambda_pow + 1.0
    #Lp3 = lambda_pow + 3.0
    #Lp4 = lambda_pow + 4.0
    #fnb = float(nbins_grains)

    #if nbins_grains == 1:
         #r_g[0] = a_0
    #else:
        #for i in range(nbins_grains):
            #r_g[i] = a_min*zeta**(-float(i+1)/fnb) * np.sqrt( Lp1/Lp3* (1.0-zeta**(Lp3/fnb))/(1.0-zeta**(Lp1/fnb)))


    #q   = np.zeros([theTable.nx[0]])
    #m   = np.zeros([theTable.nx[0]])
    #m_g = np.zeros([theTable.nx[0]])
    
    #q[:] = e
    #q[0] = -e
    #if ndims == 3:
        #qchrg = [e,-e,0.0]
        #for i in range(nion,theTable.nx[0]):
            #q[i] = qchrg[(i-7) % 3]               
    #elif ndims == 4:
        #qchrg = [0.0,e,-e]
        #for i in range(nion,theTable.nx[0]):
            #q[i] = qchrg[(i-8) % 3]               

    #m[0] = me        # e-
    #m[1] = 23.5*mp   # metallic ions
    #m[2] = 29.0*mp   # molecular ions
    #m[3] = 3.0*mp    # H3+
    #m[4] = mp        # H+
    #m[5] = 12.0*mp   # C+
    #m[6] = 4.0*mp    # He+
    #if ndims == 4:
        #m[7] = 39.098*mp # K+
        #m[8] = 22.990*mp # Na+
    #for i in range(nbins_grains):
        #m_g[i] = 4.0/3.0*np.pi*r_g[i]**3*rho_s
        #m[nion+1+3*i:nion+1+3*(i+1)] = m_g[i]
    
    
    #print q
    #print '======================='
    #print m
    
    
    ## Compute conductivities =============================
    
    ## Define magnetic field range and resolution
    #bminchimie = 1.0e-10
    #bmaxchimie = 1.0e+10
    #bchimie = 150
    ##dbchimie=(np.log10(bmaxchimie)-np.log10(bminchimie))/float(bchimie-1)
    #Barray = np.linspace(np.log10(bminchimie),np.log10(bmaxchimie),bchimie)
    
    #nchimie = theTable.nx[1]
    #tchimie = theTable.nx[2]
    #if ndims == 3:
        #xichimie = 1
        #theTable.eta_ohm  = np.zeros([nchimie,tchimie,bchimie])
        #theTable.eta_ad   = np.zeros([nchimie,tchimie,bchimie])
        #theTable.eta_hall = np.zeros([nchimie,tchimie,bchimie])
        #theTable.eta_hsig = np.zeros([nchimie,tchimie,bchimie])
    #elif ndims == 4:
        #xichimie = theTable.nx[3]
        #theTable.eta_ohm  = np.zeros([nchimie,tchimie,xichimie,bchimie])
        #theTable.eta_ad   = np.zeros([nchimie,tchimie,xichimie,bchimie])
        #theTable.eta_hall = np.zeros([nchimie,tchimie,xichimie,bchimie])
        #theTable.eta_hsig = np.zeros([nchimie,tchimie,xichimie,bchimie])

    #tau_sn      = np.zeros([theTable.nx[0]])
    #omega       = np.zeros([theTable.nx[0]])
    #sigma       = np.zeros([theTable.nx[0]])
    #phi         = np.zeros([theTable.nx[0]])
    #zetas       = np.zeros([theTable.nx[0]])
    #gamma_zeta  = np.zeros([theTable.nx[0]])
    #gamma_omega = np.zeros([theTable.nx[0]])
    #omega_bar   = np.zeros([theTable.nx[0]])
    
    #for iX in range(xichimie):
        #for iB in range(bchimie):
            #print iB,bchimie
            #for iT in range(tchimie):
                #for iH in range(nchimie):

                    #B = Barray[iB]
                    #if ndims == 3:
                        #nh = resistivite_chimie_x[0,iH,iT]  # density (.cc) of current point
                        #T  = resistivite_chimie_x[1,iH,iT]
                    #elif ndims == 4:
                        #nh = resistivite_chimie_x[0,iH,iT,iX]  # density (.cc) of current point
                        #T  = resistivite_chimie_x[1,iH,iT,iX]
                        #xi = resistivite_chimie_x[2,iH,iT,iX]
      
                    #for i in range(nion):
                        #if  i==0 : # electron
                            #if ndims == 3:
                                #sigv = 1.3e-9
                            #elif ndims == 4:
                                #sigv = 3.16e-11 * (np.sqrt(8.0*kb*1.0e-7*T/(np.pi*me*1.0e-3))*1.0e-3)**1.3
                            #tau_sn[i] = 1.0/1.16*(m[i]+2.0*mp)/(2.0*mp)*1.0/(nh/2.0*sigv)
                        
                        #else: # ions   
                        ##elif (i>=1) and (i<nion): # ions
                            #if ndims == 3:
                                #sigv = 1.69e-9
                            #elif ndims == 4:
                                #muuu=m[i]*2.0*mp/(m[i]+2.0*mp)
                                #if (i==1) or (i==2):
                                    #sigv=2.4e-9 *(np.sqrt(8.0*kb*1.0e-7*T/(np.pi*muuu*1.0e-3))*1.0e-3)**0.6
                                #elif i==3:
                                    #sigv=2.0e-9 * (np.sqrt(8.0*kb*1.0e-7*T/(np.pi*muuu*1.0e-3))*1.0e-3)**0.15
                                #elif i==4:
                                    #sigv=3.89e-9 * (np.sqrt(8.0*kb*1.0e-7*T/(np.pi*muuu*1.0e-3))*1.0e-3)**(-0.02)
                                #else:
                                    #sigv=1.69e-9
                            #tau_sn[i] = 1.0/1.14*(m[i]+2.0*mp)/(2.0*mp)*1.0/(nh/2.0*sigv)
                        
                        #omega[i] = q[i]*B/(m[i]*clight)
                        #if ndims == 3:
                            #sigma[i] = resistivite_chimie_x[i+2,iH,iT]*nh*(q[i])**2*tau_sn[i]/m[i]
                        #elif ndims == 4:
                            #sigma[i] = resistivite_chimie_x[i+3,iH,iT,iX]*nh*(q[i])**2*tau_sn[i]/m[i]
                        ##phi[i] = 0.0
                        ##zetas[i] = 0.0
                        #gamma_zeta[i] = 1.0
                        #gamma_omega[i] = 1.0
                        ##omega_bar[i] = 0.0
                        
                    #for i in range(nbins_grains):
                        
                        ## g+
                        #tau_sn[nion+1+3*i] = 1.0/1.28*(m_g[i]+2.0*mp)/(2.0*mp)*1.0/(nh/2.0*(np.pi*r_g[i]**2*(8.0*kb*T/(np.pi*2.0*mp))**0.5))
                        #omega [nion+1+3*i] = q[nion+1+3*i]*B/(m_g[i]*clight)
                        
                        ## g-
                        #tau_sn[nion+2+3*i] = tau_sn[nion+1+3*i]
                        #omega [nion+2+3*i] = q[nion+2+3*i]*B/(m_g[i]*clight)
                        
                        #if ndims == 3:
                            #sigma[nion+1+3*i] = resistivite_chimie_x[nion+3+3*i,iH,iT]*nh*(q[nion+1+3*i]**2)*tau_sn[nion+1+3*i]/m_g[i]
                            #sigma[nion+2+3*i] = resistivite_chimie_x[nion+4+3*i,iH,iT]*nh*(q[nion+2+3*i]**2)*tau_sn[nion+2+3*i]/m_g[i]
                        #elif ndims == 4:
                            #sigma[nion+1+3*i] = resistivite_chimie_x[nion+4+3*i,iH,iT,iX]*nh*(q[nion+1+3*i]**2)*tau_sn[nion+1+3*i]/m_g[i]
                            #sigma[nion+2+3*i] = resistivite_chimie_x[nion+5+3*i,iH,iT,iX]*nh*(q[nion+2+3*i]**2)*tau_sn[nion+2+3*i]/m_g[i]
                    
                    #sigP =0.0
                    #sigO =0.0
                    #sigH =0.0

                    #for i in range(theTable.nx[0]):
                        #sigP =sigP + sigma[i]
                        #sigO =sigO + sigma[i]/(1.0+(omega[i]*tau_sn[i])**2)
                        #sigH =sigH - sigma[i]*omega[i]*tau_sn[i]/(1.0+(omega[i]*tau_sn[i])**2)
                    
                    #if ndims == 3:
                        #theTable.eta_ohm [iH,iT,iB] = np.log10(1.0/sigP)                        # Ohmic
                        #theTable.eta_ad  [iH,iT,iB] = np.log10(sigO/(sigO**2+sigH**2)-1.0/sigP) # Ambipolar
                        #theTable.eta_hall[iH,iT,iB] = np.log10(abs(sigH/(sigO**2+sigH**2)))     # Hall
                        #theTable.eta_hsig[iH,iT,iB] = sigH / abs(sigH)                          # Hall sign
                    #elif ndims == 4:
                        #theTable.eta_ohm [iH,iT,iX,iB] = np.log10(1.0/sigP)                        # Ohmic
                        #theTable.eta_ad  [iH,iT,iX,iB] = np.log10(sigO/(sigO**2+sigH**2)-1.0/sigP) # Ambipolar
                        #theTable.eta_hall[iH,iT,iX,iB] = np.log10(abs(sigH/(sigO**2+sigH**2)))     # Hall
                        #theTable.eta_hsig[iH,iT,iX,iB] = sigH / abs(sigH)                          # Hall sign

    #print resistivite_chimie_x[0,:,:]
    ##theTable.dens = resistivite_chimie_x
    ##theTable.grid = (theTable.dens,theTable.tgas,theTable.trad)
