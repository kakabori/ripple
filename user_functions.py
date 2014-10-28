"""
This module implements functions that normal program users would use, 
such as plotting and exporting.
"""
from edm import *
from ripple import *
from data import *

data = Data()
rip = None
edm = ElectronDensityMap(data)

def load_data(filename):
    """Load intensity data"""
    global rip
    h, k, q, I, sigma = read_data_5_columns(filename)
    rip = M2G(h, k, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
              xM=90, A=25, f1=1.5, f2=-20, 
              rho_H1=9.91, Z_H1=20, sigma_H1=2.94,
              rho_H2=7.27, Z_H2=20, sigma_H2=1.47, 
              rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=3)
    rip.edp_par['f2'].vary = True
    rip.edp_par['rho_H1'].vary = False
    rip.edp_par['sigma_H1'].vary = False
    rip.edp_par['rho_H2'].vary = False
    rip.edp_par['sigma_H2'].vary = False
    rip.edp_par['rho_M'].vary = False
    rip.edp_par['sigma_M'].vary = False

def invert_phase(h, k):
    """Invert the sign of a phase factor for (h, k). This function
    affects EDM and EDP plots."""
    global data
    data.flip_phases(h, k)
    
def fit_edp():
    """Initiate a NLSQ fitting of the input intensity data to a model.
    After a fit is finished, the phase factors for plotting 
    will also be updated with predicted phase factors.
    """
    global rip
    global data
    rip.fit_edp()
    rip.report_edp()
    data.update(rip.h, rip.k, rip.qx, rip.qz, rip.F) 

def export_EDM(xmin=-150, xmax=150, zmin=-100, zmax=100, N=301, filename="EDM.dat"):
    """Export EDM as an ASCII file"""
    global data
    global edm
    X, Y, Z = edm.get_EDM(xmin, xmax, zmin, zmax, N, data)
    with open(filename, 'w') as f:
        f.write("x z ED\n")
        for x, y, z in zip(X, Y, Z):
            f.write("{0: 3.1f} {1: 3.1f} {2: }\n".format(x, y, z))     

def plot_EDM(xmin=-150, xmax=150, zmin=-100, zmax=100, N=301):
    """Plot electron density map."""
    # X and Y are x and z coordinates, respectively
    # Z is calclated electron densities at points (x, z)
    global data
    global edm
    X, Y, Z = edm.get_EDM(xmin, xmax, zmin, zmax, N, data)
    X.shape = (N, N)
    Y.shape = (N, N)
    Z.shape = (N, N)
    plt.figure()
    rotate_Z = ndimage.rotate(Z, 90)
    imgplot = plt.imshow(rotate_Z, extent=[xmin,xmax,zmin,zmax], cmap='gray')
    #imgplot = plt.imshow(rotate_Z, cmap='gray')
#    return imgplot 

def export_EDP_endpoints(start, end, N, filename):
    global data
    global edm
    X, Z, DIST, EDP = edm.get_EDP_endpoints(start, end, N, data)
    _export_EDP(X, Z, DIST, EDP, filename)    

def plot_EDP_endpoints(start, end, N): 
    """Plot an experimental EDP along a line connecting start 
    and end, on N points.
    """
    global data
    global edm
    X, Z, DIST, EDP = edm.get_EDP_endpoints(start, end, N, data)
    _plot_EDP(DIST, EDP) 
    
def export_EDP_angle(center, angle, length, stepsize, filename):
    global data
    global edm
    X, Z, DIST, EDP = edm.get_EDP_angle(center, angle, length, stepsize, data)
    _export_EDP(X, Z, DIST, EDP, filename)
     
def plot_EDP_angle(center, angle, length, stepsize):
    global data
    global edm
    X, Z, DIST, EDP = edm.get_EDP_angle(center, angle, length, stepsize, data)
    _plot_EDP(DIST, EDP) 
    
def _export_EDP(X, Z, DIST, EDP, filename):
    with open(filename, 'w') as f:
        f.write("x z dist ED\n")
        for x, z, dist, edp in zip(X, Z, DIST, EDP):
            f.write("{0: 3.1f} {1: 3.1f} {2: 3.1f} {3: }\n".format(x, z, dist, edp))

def _plot_EDP(DIST, EDP):
    plt.figure()
    plt.plot(DIST, EDP)
