import osiris as pp
import sys

fname = ""
if len(sys.argv) > 2:
    fname = sys.argv[2]
if len(sys.argv) > 1:
    nout = int(sys.argv[1])
else:
    nout = -1

# Load data
mydata = pp.RamsesData(nout=nout,center="max:density",scale="au")

# Density vs T
if(len(fname) > 0):
    mydata.plot_histogram("log_rho","log_B",fname=fname)
else:
    mydata.plot_histogram("log_rho","log_B",block=True)
