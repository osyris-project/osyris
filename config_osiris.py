#===================================================================================
# Define default values so that you don't have to specify them every time.
#===================================================================================
default_values = {
    "nout"   : 1    ,
    "lmax"   : 0    ,
    "center" : None ,
    "dx"     : 0.0  ,
    "dy"     : 0.0  ,
    "dz"     : 0.0  ,
    "scale"  : "au" ,
    "verbose": False,
    "path"   : "",
    }

#=======================================================================================
# Common variables
#=======================================================================================
constants = {"cm"  : 1.0         ,\
             "au"  : 1.495980e+13,\
             "pc"  : 3.085678e+18,\
             "yr"  : 365.25*86400.0,\
             "kyr" : 365.25*86400.0*1000.0,\
             "msun": 1.9889e33}

#===================================================================================
# Here are some additional variables that are to be computed every time data is
# loaded.
#===================================================================================
def additional_variables(holder):
    
    # Velocity field (in case conservative variables are dumped)
    holder.new_field(name="velocity_x",operation="momentum_x/density",unit="cm/s",label="velocity_x")
    holder.new_field(name="velocity_y",operation="momentum_y/density",unit="cm/s",label="velocity_y")
    holder.new_field(name="velocity_z",operation="momentum_z/density",unit="cm/s",label="velocity_z")
    
    # Magnetic field
    holder.new_field(name="B_x",operation="0.5*(B_left_x+B_right_x)",unit="G",label="B_x")
    holder.new_field(name="B_y",operation="0.5*(B_left_y+B_right_y)",unit="G",label="B_x")
    holder.new_field(name="B_z",operation="0.5*(B_left_z+B_right_z)",unit="G",label="B_x")
    holder.new_field(name="B",operation="np.sqrt(B_x**2+B_y**2+B_z**2)",unit="G",label="B")
    
    # Mass and radius
    holder.new_field(name="r",operation="np.sqrt(x**2 + y**2 + z**2)",unit="cm",label="Radius")
    holder.new_field(name="mass",operation="density*((dx*constants[\""+holder.info["scale"]+"\"])**3)/constants[\"msun\"]",unit="Msun",label="Mass")
    
    # Commonly used log quantities
    holder.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    holder.new_field(name="log_T",operation="np.log10(temperature)",unit="K",label="log(T)")
    holder.new_field(name="log_B",operation="np.log10(B)",unit="G",label="log(B)")
    holder.new_field(name="log_r",operation="np.log10(r)",unit="cm",label="log(Radius)")
    holder.new_field(name="log_m",operation="np.log10(mass)",unit="g",label="log(Mass)")
    
    #========================== ADD YOUR VARIABLES HERE ============================

    return