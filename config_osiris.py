#===================================================================================
# Define default values so that you don't have to specify them every time.
#===================================================================================
default_values = {
    "nout"      : 1    ,
    "lmax"      : 0    ,
    "center"    : None ,
    "dx"        : 0.0  ,
    "dy"        : 0.0  ,
    "dz"        : 0.0  ,
    "scale"     : "au" ,
    "verbose"   : False,
    "path"      : ""   ,
    "variables" : []   ,
    "var_names" : ["density","velocity_x","velocity_y","velocity_z","B_left_x",\
                   "B_left_y","B_left_z","B_right_x","B_right_y","B_right_z",\
                   "thermal_pressure","radiative_energy_1","passive_scalar_1",\
                   "passive_scalar_2","passive_scalar_3","passive_scalar_4",\
                   "temperature"]
}

#=======================================================================================
# Common variables
#=======================================================================================
constants = {
    "cm"  : 1.0,
    "au"  : 1.495980e+13,
    "pc"  : 3.085678e+18,
    "yr"  : 365.25*86400.0,
    "kyr" : 365.25*86400.0*1000.0,
    "msun": 1.9889e33
}

#===================================================================================
# Here are some additional variables that are to be computed every time data is
# loaded.
#===================================================================================
def additional_variables(holder):
    
    # Velocity field (in case conservative variables are dumped)
    holder.new_field(name="velocity_x",operation="momentum_x/density",unit="cm/s",label="velocity_x",verbose=False)
    holder.new_field(name="velocity_y",operation="momentum_y/density",unit="cm/s",label="velocity_y",verbose=False)
    holder.new_field(name="velocity_z",operation="momentum_z/density",unit="cm/s",label="velocity_z",verbose=False)
    
    # Magnetic field
    holder.new_field(name="B_x",operation="0.5*(B_left_x+B_right_x)",unit="G",label="B_x",verbose=False)
    holder.new_field(name="B_y",operation="0.5*(B_left_y+B_right_y)",unit="G",label="B_x",verbose=False)
    holder.new_field(name="B_z",operation="0.5*(B_left_z+B_right_z)",unit="G",label="B_x",verbose=False)
    holder.new_field(name="B",operation="np.sqrt(B_x**2+B_y**2+B_z**2)",unit="G",label="B",verbose=False)
    
    # Mass and radius
    holder.new_field(name="r",operation="np.sqrt(x**2 + y**2 + z**2)",unit=holder.info["scale"],label="Radius",verbose=False)
    holder.new_field(name="mass",operation="density*((dx*"+str(constants[holder.info["scale"]])+")**3)/"+str(constants["msun"]),unit="Msun",label="Mass",verbose=False)
    
    # Commonly used log quantities
    holder.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)",verbose=False)
    holder.new_field(name="log_T",operation="np.log10(temperature)",unit="K",label="log(T)",verbose=False)
    holder.new_field(name="log_B",operation="np.log10(B)",unit="G",label="log(B)",verbose=False)
    holder.new_field(name="log_r",operation="np.log10(r)",unit=holder.info["scale"],label="log(Radius)",verbose=False)
    holder.new_field(name="log_m",operation="np.log10(mass)",unit="g",label="log(Mass)",verbose=False)
    
    #========================== ADD YOUR VARIABLES HERE ============================

    return
