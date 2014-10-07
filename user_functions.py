"""
This module implements functions that normal program users would use, 
such as plotting and exporting.
"""
import edm as *

data = Data()
rip = None

def load_data():
    global rip
    
    
def edp_fit():
    global rip
    rip.edp_fit()
    rip.report_edp()

def export_EDM(xmin, xmax, zmin, zmax, N, data, filename):
    """Export EDM as an ASCII file"""
    edm = ElectronDensityMap(data)
    X, Y, Z = edm.get_EDM(xmin, xmax, zmin, zmax, N)
    with open(filename, 'w') as f:
        f.write("x z ED\n")
        for x, y, z in zip(X, Y, Z):
            f.write("{0: 3.1f} {1: 3.1f} {2: }\n".format(x, y, z))     

def plot_EDM(xmin, xmax, zmin, zmax, N, data):
    """Plot electron density map."""
    # X and Y are x and z coordinates, respectively
    # Z is calclated electron densities at points (x, z)
    edm = ElectronDensityMap(data)
    X, Y, Z = edm.get_EDM(xmin, xmax, zmin, zmax, N)
    X.shape = (N, N)
    Y.shape = (N, N)
    Z.shape = (N, N)
    plt.figure()
    rotate_Z = ndimage.rotate(Z, 90)
    imgplot = plt.imshow(rotate_Z, extent=[-150,150,-100,100], cmap='gray')
#   return imgplot 

def export_EDP_endpoints(start, end, N, data, filename):
    edm = ElectronDensityMap(data)
    X, Z, DIST, EDP = edm.get_EDP_endpoints(start, end, N)
    _export_EDP(X, Z, DIST, EDP, filename)    

def plot_EDP_endpoints(start, end, N, data): 
    """Plot an experimental EDP along a line connecting start 
    and end, on N points.
    """
    edm = ElectronDensityMap(data)
    X, Z, DIST, EDP = edm.get_EDP_endpoints(start, end, N):
    _plot_EDP(DIST, EDP) 
    
def export_EDP_angle(center, angle, length, stepsize, data, filename):
    edm = ElectronDensityMap(data)
    X, Z, DIST, EDP = edm.get_EDP_angle(center, angle, length, stepsize, data)
    _export_EDP(X, Z, DIST, EDP, filename)
     
def plot_EDP_angle(center, angle, length, stepsize, data):
    edm = ElectronDensityMap(data)
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
