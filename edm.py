import numpy as np
import matplotlib as ppl
import matplotlib.pyplot as plt

class ElectronDensityMap(object):
    """Implements electron density map (EDM) related methods."""      
    def __init__(self, data=None):
        """Input is Data object"""
        if data is not None:
            self.data = data
        
    def update_data(self, data):
        """Input is Data object"""
        self.data = data
                   
    def plot_EDM(self, xmin, xmax, zmin, zmax, N, data):
        """Plot electron density map."""
        # X and Y are x and z coordinates, respectively
        # Z is calclated electron densities at points (x, z)
        X, Y, Z = self.get_EDM(xmin, xmax, zmin, zmax, N, data)
        X.shape = (N, N)
        Y.shape = (N, N)
        Z.shape = (N, N)
        plt.figure()
        rotate_Z = ndimage.rotate(Z, 90)
        imgplot = plt.imshow(rotate_Z, extent=[-150,150,-100,100], cmap='gray')
#        return imgplot            

    
    def export_EDM(self, xmin, xmax, zmin, zmax, N, filename, data):
        """Export EDM as an ASCII file"""
        X, Y, Z = self.get_EDM(xmin, xmax, zmin, zmax, N)
        with open(filename, 'w') as f:
            f.write("x z ED\n")
            for x, y, z in zip(X, Y, Z):
                f.write("{0: 3.1f} {1: 3.1f} {2: }\n".format(x, y, z))        
                    
    def get_EDM(self, xmin, xmax, zmin, zmax, N, data):
        """Fourier-reconstruct a 2D map of the electron density profile and return
        as Z on (X,Y) grids. Calculate EDP at N points along x and N points along z. 
        The units are in Angstrom.
    
        output: X, Y, Z, each being numpy array
        """
        rho_xz = []
        xgrid = np.linspace(xmin, xmax, num=N)
        zgrid = np.linspace(zmin, zmax, num=N)
        F = data.form_factors()
        qx, qz = data.qx_qz()
        for x in xgrid:
            for z in zgrid:
                tmp = F * np.cos(qx*x+qz*z)
                rho_xz.append([x, z, tmp.sum(axis=0)])
        rho_xz = np.array(rho_xz, float)  
        X, Y, Z= rho_xz[:,0], rho_xz[:,1], rho_xz[:,2]
        return X, Y, Z
    
    def get_EDP_endpoints(self, start, end, N, data):
        x0, z0 = start
        x1, z1 = end
        xpoints = np.linspace(x0, x1, N)
        zpoints = np.linspace(z0, z1, N)
        X, Z, DIST, EDP = self._calc_EDP(xpoints, zpoints, start)
        return X, Z, DIST, EDP
        
    
    def plot_EDP_endpoints(self, start, end, N, F, filename=None, data): 
        """Plot an experimental EDP along a line connecting start 
        and end, on N points. If filename is specified, export an
        ASCII file instead.
        """
        x0, z0 = start
        x1, z1 = end
        xpoints = np.linspace(x0, x1, N)
        zpoints = np.linspace(z0, z1, N)
        return self._plot_EDP(xpoints, zpoints, start, filename)
    
    def get_EDP_angle(self, center, angle, length, stepsize, data):
        x, z = center
        N = length/stepsize + 1
        angle = angle*pi/180 
        if angle==0:
            # If angle is zero, the slope is infinite. 
            # In this case, x is constant.
            xpoints = x * np.ones(N)
            zpoints = np.linspace(z-length/2, z+length/2, N)
        else:
            slope = 1 / tan(angle)
            intercept = z - slope*x
            xpoints = np.linspace(x-length*sin(angle)/2, x+length*sin(angle)/2, N)
            zpoints = slope * xpoints + intercept
        X, Z, DIST, EDP = self._calc_EDP(xpoints, zpoints, center)
        return X, Z, DIST, EDP          
            
    def plot_EDP_angle(self, center, angle, length, stepsize, filename=None, data):
        x, z = center
        N = length/stepsize + 1
        angle = angle*pi/180 
        if angle==0:
            # If angle is zero, the slope is infinite. 
            # In this case, x is constant.
            xpoints = x * np.ones(N)
            zpoints = np.linspace(z-length/2, z+length/2, N)
        else:
            slope = 1 / tan(angle)
            intercept = z - slope*x
            xpoints = np.linspace(x-length*sin(angle)/2, x+length*sin(angle)/2, N)
            zpoints = slope * xpoints + intercept    
        return self._plot_EDP(xpoints, zpoints, center, filename)     
    
    def _plot_EDP(self, xarray, zarray, center, filename, data):
        F = self.data.form_factors()
        X, Z, DIST, EDP = self._calc_EDP(xarray, zarray, center)
        if filename is None:
            plt.figure()
            plt.plot(DIST, EDP)
        else:
            with open(filename, 'w') as f:
                f.write("x z dist ED\n")
                for x, z, dist, edp in zip(X, Z, DIST, EDP):
                    f.write("{0: 3.1f} {1: 3.1f} {2: 3.1f} {3: }\n".format(x, z, dist, edp))
        return X, Z, DIST, EDP
        
    def _calc_EDP(self, xpoints, zpoints, center, data):
        xM, z0 = center
        rho = []
        F = self.data.form_factors()
        qx, qz = self.data.qx_qz()
        for x, z in zip(xpoints, zpoints):
            tmp = F * np.cos(qx*x+qz*z)
            dist = np.sign(z-z0)*np.sqrt((x-xM)**2 + (z-z0)**2)
            rho.append([x, z, dist, tmp.sum(axis=0)])
        rho = np.array(rho, float)
        X, Z, DIST, EDP = rho[:,0], rho[:,1], rho[:,2], rho[:,3]   
        return X, Z, DIST, EDP     
  
    def export_headgroup_positions(self, lambda_r, D, A, xM, filename):
        stepsize = 1
        xmin, xmax = -lambda_r, lambda_r
        length = xmax - xmin
        x_array = np.linspace(xmin, xmax, int(math.ceil(length)/stepsize+1))
        z_low_list, ed_low_list, z_up_list, ed_up_list = [], [], [], []
        for x in x_array:
            z_low, ed_low, z_up, ed_up = self.find_headgroup(x, lambda_r, D, A, xM)
            z_low_list.append(z_low)
            ed_low_list.append(ed_low)
            z_up_list.append(z_up)
            ed_up_list.append(ed_up)
        with open(filename, 'w') as f:
            f.write("x z_lower ED_lower z_upper ED_upper\n")
            for a, b, c, d, e in zip(x_array, z_low_list, ed_low_list, z_up_list, ed_up_list):
                f.write("{0: 3.1f} {1: 3.1f} {2: 6.1f} {3: 3.1f} {4: 6.1f}\n".format(a, b, c, d, e))
      
    def export_methyl_positions(self, lambda_r, D, A, xM, filename):
        stepsize = 1
        xmin, xmax = -lambda_r, lambda_r
        length = xmax - xmin
        x_array = np.linspace(xmin, xmax, int(math.ceil(length)/stepsize+1))
        z_list, ed_list = [], []
        for x in x_array:
            z, ed = self.find_methyl(x, lambda_r, D, A, xM)
            z_list.append(z)
            ed_list.append(ed)
        with open(filename, 'w') as f:
            f.write("x z ED\n")
            for a, b, c in zip(x_array, z_list, ed_list):
                f.write("{0: 3.1f} {1: 3.1f} {2: 6.1f}\n".format(a, b, c))    
    
    def find_headgroup(self, x, lambda_r, D, A, xM):
        """Return the z position of maximum electron density along a vertical line
        at x and corresponding electron density. For a normal EDP, this should 
        correspond to the headgroup position.
        """
        z0 = where_in_sawtooth(x, lambda_r, A, xM)
        # D*10+1 is the number of points, every 0.1 Angstrom
        z_lower = np.linspace(z0-D/2, z0, D/2*10+1)
        z_upper = np.linspace(z0, z0+D/2, D/2*10+1) 
        x_array = np.zeros(D/2*10+1) + x
        edp_lower = self._get_electron_density(x_array, z_lower)
        edp_upper = self._get_electron_density(x_array, z_upper)
        return (z_lower[np.argmax(edp_lower)], np.amax(edp_lower), 
                z_upper[np.argmax(edp_upper)], np.amax(edp_upper))
  
    def find_methyl(self, x, lambda_r, D, A, xM):
        """Return the z position of minimum electron density along a vertical line
        at x and corresponding electron density. For a normal EDP, this should
        correspond to the terminal methyl group position.
        """
        z0 = where_in_sawtooth(x, lambda_r, A, xM)
        # +/-10 Angstrom from the sawtooth should be enough to find the mid-plane
        z_array = np.linspace(z0-10, z0+10, 201)
        x_array = np.zeros(201) + x
        edp = self._get_electron_density(x_array, z_array)
        return (z_array[np.argmin(edp)], np.amin(edp))
  
    def _get_electron_density(self, x_array, z_array):
        """Return electron density calculated at points specified by 
        x_array and z_array
        """
        if x_array.size != z_array.size:
            print("length of x must be equal to length of z")
            return
        tmp = np.zeros(x_array.size)
        for F, qx, qz in zip(self.F, self.qx, self.qz):
            tmp = tmp + F * cos(qx*x_array+qz*z_array)
        return tmp
        

###############################################################################
