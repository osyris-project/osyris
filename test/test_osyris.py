import pytest
import osyris
import numpy as np
import matplotlib.pyplot as plt


# @pytest.fixture(scope="session", autouse=True)
# def load_ramses_output(request):
#        mydata = osyris.RamsesData(nout=71, center="max:density",
#                               scale="au")
#     # prepare something ahead of all tests
#     request.addfinalizer(finalizer_function)


def test_simple_2d_slice():

    mydata = osyris.RamsesData(71, scale="au")
    osyris.plot_slice(mydata.log_rho, direction="z",
                      vec=mydata.velocity, dx=100)


def test_2d_histogram():

    mydata = osyris.RamsesData(71)
    osyris.plot_histogram(mydata.log_rho, mydata.log_T,
                          scalar_args={"cmap": "log"})


def test_new_variable_field_and_scatter_plot():

    mydata = osyris.RamsesData(nout=71, center=[0.5, 0.5, 0.5], scale="au")
    mydata.new_field(name="log_vel",
                     operation="np.log10(np.sqrt(velocity_x**2 +"
                               "velocity_y**2 + velocity_z**2))",
                     unit="cm/s", label="log(Velocity)")
    osyris.plot_histogram(mydata.log_rho, mydata.log_vel, scatter=mydata.log_T,
                          outline=True,
                          scatter_args={"iskip": 100, "cmap": "gnuplot"})


def test_center_max_density_slice_with_streamlines():

    osyris.config.default_values["time_unit"] = "yr"
    mydata = osyris.RamsesData(nout=71, center="max:density",
                               scale="au")
    osyris.plot_slice(mydata.log_rho, direction="yxz",
                      stream=mydata.B, dx=100,
                      stream_args={"cmap": "log,jet"})


def test_arbitrary_angle_slice_with_coloured_velocity_vectors():

    mydata = osyris.RamsesData(nout=71, center="max:density",
                               scale="au")
    mydata.new_field(name="vkms", operation="velocity/1.0e5",
                     unit="km/s", label="Velocity")
    osyris.plot_slice(mydata.log_rho, direction=[-1, 1, 1],
                      vec=mydata.vkms, dx=100,
                      vec_args={"cmap": "jet", "vskip": 4, "cbar": True})


def test_automatic_slice_orientation_with_ang_mom():

    mydata = osyris.RamsesData(nout=71, center="max:density",
                               scale="au")
    osyris.plot_slice(mydata.log_rho, direction="auto:top",
                      vec=mydata.velocity, dx=100)
    osyris.plot_slice(mydata.log_rho, direction="auto:side",
                      vec=mydata.velocity, dx=100)


def test_subplots_slices_and_contours():

    # Create figure
    fig = plt.figure(figsize=(15, 5.25))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Load data
    mydata = osyris.RamsesData(nout=71, center="av:density>1.0e-10",
                               scale="au", dx=20.0, dy=20.0, dz=20.0)

    # Create velocity field in km/s
    mydata.new_field(name="vkms", operation="velocity/1.0e5",
                     unit="km/s", label="Velocity")

    # Define region to plot
    dx = 15.0

    # Left plot: coloured density slice with overlayed contours
    osyris.plot_slice(mydata.log_rho, direction="z", dx=dx,
                      axes=ax1,
                      scalar_args={"extend": "both",
                                   "vmin": -14.0,
                                   "vmax": -9.0,
                                   "nc": 40},
                      vec=mydata.vkms,
                      vec_args={"vscale": 2.0,
                                "vkey_pos": [0.65, 0.1]},
                      contour=mydata.log_rho,
                      contour_args={"levels": [-12.0, -11.0, -9.0],
                                    "colors": ("yellow", "k", "lime"),
                                    "linewidths": [2, 5, 2],
                                    "linestyles": ["solid", "dashed", "solid"],
                                    "cmap": None,
                                    "cbar": False},
                      title="My title")

    # Right plot: temperature slice with AMR levels
    osyris.plot_slice(mydata.log_T, direction="z", dx=dx,
                      axes=ax2, title="",
                      scalar_args={"cmap": "hot"},
                      contour=mydata.level,
                      contour_args={"fmt": "%i",
                                    "colors": "lightgray",
                                    "cmap": None,
                                    "levels": range(12, 20),
                                    "label": True,
                                    "cbar": False})

    # Save figure to pdf file
    fig.savefig("demo008.png")


def test_plot_subset_of_disk_cells():

    # Load data
    mydata = osyris.RamsesData(nout=71, center="max:density",
                               scale="au", dx=100.0, dy=100.0,
                               dz=100.0)

    mydata.new_field(name="log_rho_disk",
                     values=np.where(np.logical_and(
                         mydata.get("log_rho", only_leafs=False) > -12.5,
                         mydata.get("log_rho", only_leafs=False) < -11.0),
                         mydata.get("log_rho", only_leafs=False), np.NaN),
                     label="Disk density")

    osyris.plot_slice(mydata.log_rho_disk, direction="z", dx=50)

    # Now print disk mass: 2 different ways
    # Method 1:
    cube = np.where(np.logical_and(
                   mydata.get("log_rho") > -12.5,
                   mydata.get("log_rho") < -11.0))
    mcore1 = np.sum(mydata.get("mass")[cube])
    # Method 2:
    mydata.new_field(name="disk_mass",
                     values=np.where(np.logical_and(
                         mydata.get("log_rho", only_leafs=False) > -12.5,
                         mydata.get("log_rho", only_leafs=False) < -11.0),
                         mydata.get("mass", only_leafs=False), np.NaN),
                     label="Disk mass")
    mcore2 = np.nansum(mydata.get("disk_mass"))
    print("Disk mass: %.3e Msun ; %.3e Msun"%(mcore1, mcore2))


def test_difference_between_two_snapshots():

    # Read data from 2 snapshots
    mydata1 = osyris.RamsesData(71, scale="au")
    mydata2 = osyris.RamsesData(201, scale="au")

    # Extract log(density) slices by copying data into structures
    slice1 = osyris.plot_slice(mydata1.log_rho, direction="z",
                               dx=100, plot=False, copy=True)
    slice2 = osyris.plot_slice(mydata2.log_rho, direction="z",
                               dx=100, plot=False, copy=True)

    # Get coordinates
    x = slice1[0]
    y = slice1[1]

    # Get densities
    rho1 = slice1[2]
    rho2 = slice2[2]

    # Density difference
    diff = (rho1-rho2)/rho2

    # Create figure
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    im1 = ax1.contourf(x, y, diff, cmap='RdBu',
                       levels=np.linspace(-0.12, 0.12, 31))
    ax1.set_aspect("equal")
    cb = plt.colorbar(im1, ax=ax1)
    cb.ax.set_ylabel("Relative difference")

    fig.savefig("diff.png", bbox_inches="tight")


@pytest.mark.skip(reason="ISM tables are in a separate repository")
def test_ism_tables():

    # Create figure
    fig = plt.figure(figsize=(20, 5))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Load data
    mydata = osyris.RamsesData(nout=71, scale="au",
                               verbose=True)

    # Create a properly dimensioned field for internal energy
    conv = (mydata.info["unit_l"] / mydata.info["unit_t"])**2
    mydata.new_field(name="internal_energy",
                     operation="passive_scalar_4*density*" + str(conv),
                     unit="erg/cm3",
                     label="Internal energy")

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
    ax3.text(-16.0, 0.0, "Ambipolar", va="center", ha="center")
    ax3.text(-12.0, -4.0, "Ohmic", va="center", ha="center")

    fig.savefig("ism_tables.png", bbox_inches="tight")


def test_make_a_thick_slice():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    osyris.plot_column_density(scalar=mydata.density,
                               direction="z", vec=mydata.velocity,
                               dx=100, dz=20,
                               scalar_args={"cmap": "log"},
                               nz=5, summed=False)


def test_slice_above_the_origin():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    osyris.plot_slice(scalar=mydata.density, direction="z",
                      vec=mydata.velocity, dx=100, origin=[0, 0, 5],
                      scalar_args={"cmap": "log"}, fname="demo013.png")


def test_histogram_with_mass_colormap():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    osyris.plot_histogram(mydata.log_rho, mydata.log_T, scalar=mydata.mass,
                          summed=True, scalar_args={"cmap": "magma_r,log"},
                          outline=True, fname="demo014.png")


def test_demo_with_subplots():

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
                                        "levels": range(5, 20), "cbar": False})

    # Create new field with log of velocity
    mydata.new_field(name="log_vel",
                     operation="np.log10(np.sqrt("
                               "velocity_x**2+velocity_y**2+velocity_z**2))",
                     unit="cm/s", label="log(Velocity)")

    # Density vs log_vel in scatter mode with a grey outline
    osyris.plot_histogram(mydata.log_rho, mydata.log_vel, axes=ax2,
                          scatter=mydata.log_T, outline=True,
                          scatter_args={"iskip": 100, "cmap": "gnuplot"})

    # x,z density slice with B field streamlines
    osyris.plot_slice(mydata.density, direction="yxz", stream=mydata.B,
                      dx=100, axes=ax3, scalar_args={"cmap": "log"})
    # x,y density slice with velocity vectors in color
    osyris.plot_slice(scalar=mydata.log_rho, direction="z",
                      vec=mydata.velocity, dx=100, axes=ax4,
                      vec_args={"cmap": "seismic", "vskip": 4})
    # x,y temperature slice with velocity vectors
    osyris.plot_slice(mydata.log_T, direction="z", vec=mydata.velocity,
                      dx=100, axes=ax5, scalar_args={"cmap": "hot"},
                      contour=mydata.level,
                      contour_args={"fmt": "%i", "label": True,
                                    "colors": "w", "cmap": None,
                                    "levels": range(9, 17)})

    # Now update values with later snapshot
    mydata.update_values(201)
    # Re-plot x,y density slice with velocity vectors
    osyris.plot_slice(mydata.log_rho, direction="auto:top",
                      vec=mydata.velocity, dx=100, axes=ax6)

    fig.savefig("demo015.png", bbox_inches="tight")


def test_color_slice_vectors_with_custom_field():

    mydata = osyris.RamsesData(nout=71, center="max:density", scale="au")
    mydata.new_field(name="vkms", operation="velocity/1.0e5",
                     unit="km/s", label="Velocity")
    osyris.plot_slice(scalar=mydata.log_rho, direction="yxz",
                      vec=mydata.B, dx=100,
                      scalar_args={"cmap": "Blues"},
                      vec_args={"cmap": "YlOrRd", "colors": mydata.vkms,
                                "normalize_arrows": True, "vkey": False,
                                "scale": 25.0, "cbar": True, "width": 0.01,
                                "headwidth": 1, "headlength": 0},
                      fname="demo016.png")


def test_radial_profile():

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
    re = np.linspace(rmin, rmax, nr+1)
    log_r = np.zeros([nr])
    for i in range(nr):
        log_r[i] = 0.5*(re[i]+re[i+1])

    # Modify r values so that the central cell is not "-inf"
    r = np.where(np.isinf(mydata.log_r.values), -2.0, mydata.log_r.values)

    # Bin the data in radial bins
    z0, edges = np.histogram(r, bins=re)
    z1, edges = np.histogram(r, bins=re, weights=mydata.density.values)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_mean = np.log10(z1 / z0)

    # Overlay profile
    ax.plot(log_r, rho_mean, color="r", lw=3, label="Mean profile")
    ax.legend()

    fig.savefig("demo017.png", bbox_inches='tight')
