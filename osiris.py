import config_osiris as conf
import load_ramses_data
import numpy as np

#=======================================================================================
np.seterr(divide="ignore",invalid="ignore") # Ignore divide by zero warnings
#=======================================================================================

# Define one class per type of snapshot to read ============================================

# Ramses data format =======================================================================
class RamsesData(load_ramses_data.LoadRamsesData):
 
    def __init__(self,nout=conf.default_values["nout"],lmax=conf.default_values["lmax"],\
                 center=conf.default_values["center"],dx=conf.default_values["dx"],\
                 dy=conf.default_values["dy"],dz=conf.default_values["dz"],\
                 scale=conf.default_values["scale"],verbose=conf.default_values["verbose"],\
                 nmaxcells=conf.default_values["nmaxcells"],path=conf.default_values["path"],\
                 variables=conf.default_values["variables"]):
        
        load_ramses_data.LoadRamsesData.__init__(self,nout=nout,lmax=lmax,center=center,\
                     dx=dx,dy=dy,dz=dz,scale=scale,verbose=verbose,nmaxcells=nmaxcells,\
                     path=path,variables=variables)
        
        return

# Heracles data format =====================================================================
