import osyris
import numpy as np
import matplotlib.pyplot as plt

def test_demo011():

    # Create figure
    fig = plt.figure(figsize=(20, 5))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Load data
    mydata = osyris.RamsesData(nout=71, scale="au",
                               verbose=True)

    conv = (mydata.info["unit_l"] / mydata.info["unit_t"])**2

    mydata.new_field(name="internal_energy",
                     operation="passive_scalar_4*density*" + str(conv),
                     unit="erg/cm3",
                     label="Internal energy")

    mydata.print_info()

    # Read EOS table and plot sound speed histogram
    osyris.ism_physics.get_eos(mydata, fname="SCvH_eos.dat")
    mydata.new_field(name="log_cs", operation="np.log10(cs_eos)",
                     label="log(cs)", unit="cm/s")
    osyris.plot_histogram(mydata.log_rho, mydata.log_cs,
                          scalar_args={"cmap": "log,Greens",
                                       "cbar": False},
                          outline=True, axes=ax1,
                          title="Equation of state")

    # Read opacity table and plot Rosseland mean opacity
    osyris.ism_physics.get_opacities(mydata)
    mydata.new_field(name="log_kr", operation="np.log10(kappa_r)",
                     label="log(Kr)", unit="cm2/g")
    osyris.plot_histogram(mydata.log_T, mydata.log_kr,
                          scalar_args={"cmap": "log,Blues",
                                       "cbar": False},
                          outline=True, axes=ax2,
                          title="Opacities")

    # Read resistivity table and plot Ohmic and Ambipolar
    osyris.ism_physics.get_resistivities(mydata)
    mydata.new_field(name="log_etaO", operation="np.log10(eta_ohm)",
                     label="log(etaO)")
    mydata.new_field(name="log_etaA", operation="np.log10(eta_ad)",
                     label="log(etaA)")
    osyris.plot_histogram(mydata.log_rho, mydata.log_etaO,
                          scalar_args={"cmap": "log,Greys",
                                       "cbar": False},
                          outline=True, axes=ax3, title="")
    osyris.plot_histogram(mydata.log_rho, mydata.log_etaA,
                          scalar_args={"cmap": "log,Reds",
                                       "cbar": False},
                          outline=True, axes=ax3,
                          title="Resistivities")
    ax3.set_ylabel("log(eta) [s]")
    ax3.text(-16.0,0.0, "Ambipolar", va="center", ha="center")
    ax3.text(-12.0,-4.0, "Ohmic", va="center", ha="center")

    fig.savefig("ism_tables.png", bbox_inches="tight")


def test_demo012():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    osyris.plot_column_density(scalar=mydata.density,
                               direction="z", vec=mydata.velocity,
                               dx=100, dz=20,
                               scalar_args={"cmap": "log"},
                               nz=5, summed=False, fname="demo012.png")


def test_demo013():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    osyris.plot_slice(scalar=mydata.density, direction="z",
                      vec=mydata.velocity, dx=100, origin=[0,0,5],
                      scalar_args={"cmap": "log"}, fname="demo013.png")


def test_demo014():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    osyris.plot_histogram(mydata.log_rho, mydata.log_T, mydata.mass,
                          summed=True, scalar_args={"cmap": "magma_r,log"},
                          outline=True, fname="demo014.png")


def test_demo015():

    # Load data
    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")

    # Create figure
    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    # Density vs B field with AMR level contours
    osyris.plot_histogram(mydata.log_rho, mydata.log_B, axes=ax1,
                          scalar=True, scalar_args={"cmap": "log,YlGnBu"},
                          contour=mydata.level,
                          contour_args={"fmt": "%i", "label": True,
                                        "colors": "k", "cmap": None,
                                        "levels": range(5,20), "cbar": False})

    # Create new field with log of velocity
    mydata.new_field(name="log_vel",
                     operation="np.log10(np.sqrt("
                               "velocity_x**2+velocity_y**2+velocity_z**2))",
                     unit="cm/s", label="log(Velocity)")

    # Density vs log_vel in scatter mode with a grey outline
    osyris.plot_histogram(mydata.log_rho, mydata.log_vel, axes=ax2,
                          scatter=mydata.log_T, outline=True,
                          scatter_args={"iskip": 100, "cmap": "gnuplot"})

    #x,z density slice with B field streamlines
    osyris.plot_slice(mydata.density, direction="yxz", stream=mydata.B,
                      dx=100, axes=ax3, scalar_args={"cmap": "log"})
    # x,y density slice with velocity vectors in color
    osyris.plot_slice(scalar=mydata.log_rho, direction="z",
                      vec=mydata.velocity, dx=100, axes=ax4,
                      vec_args={"cmap": "seismic", "vskip": 4})
    # x,y temperature slice with velocity vectors
    osyris.plot_slice(mydata.log_T, direction="z", vec=mydata.velocity,
                      dx=100, axes=ax5, scalar_args={"cmap":"hot"},
                      contour=mydata.level,
                      contour_args={"fmt": "%i", "label": True,
                                    "colors": "w", "cmap": None,
                                    "levels": range(9,17)})

    # Now update values with later snapshot
    mydata.update_values(201)
    # Re-plot x,y density slice with velocity vectors
    osyris.plot_slice(mydata.log_rho, direction="auto:top",
                      vec=mydata.velocity, dx=100, axes=ax6)

    fig.savefig("demo015.png", bbox_inches="tight")


def test_demo016():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    mydata.new_field(name="vkms", operation="velocity/1.0e5",
                     unit="km/s", label="Velocity")
    osyris.plot_slice(scalar=mydata.log_rho, direction="yxz",
                      vec=mydata.B, dx=100,
                      scalar_args={"cmap": "Blues"},
                      vec_args={"cmap":"YlOrRd", "colors": mydata.vkms,
                                "normalize_arrows": True, "vkey": False,
                                "scale": 25.0, "cbar": True, "width": 0.01,
                                "headwidth": 1, "headlength":0},
                      fname="demo016.png")


def test_demo017():

    # Load data
    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")

    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Make scatter plot as radial profile
    osyris.plot_histogram(mydata.log_r, mydata.log_rho, scatter=True,
                          scatter_args={"iskip": 100, "c": "grey"},
                          axes=ax)

    # Now overlay mean profile

    # Define min and max range
    rmin = -1.0
    rmax = 4.0

    # Number of points
    nr = 200

    # Radial bin edges and centers
    re = np.linspace(rmin,rmax,nr+1)
    log_r = np.zeros([nr])
    for i in range(nr):
        log_r[i] = 0.5*(re[i]+re[i+1])

    # Modify r values so that the central cell is not "-inf"
    r = np.where(np.isinf(mydata.log_r.values),-2.0,mydata.log_r.values)

    # Bin the data in radial bins
    z0, edges = np.histogram(r, bins=re)
    z1, edges = np.histogram(r, bins=re, weights=mydata.density.values)
    rho_mean = np.log10(z1 / z0)

    # Overlay profile
    ax.plot(log_r, rho_mean, color="r", lw=3, label="Mean profile")
    ax.legend()

    fig.savefig("demo017.png", bbox_inches='tight')


if __name__ == '__main__':
    test_demo017()
