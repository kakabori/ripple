import sys
import numpy as np
from numpy import pi, sin, cos, tan, exp, sqrt
import matplotlib as ppl
import scipy.optimize
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit
from scipy import ndimage
import math

# A module-level global variable
wavelength = 1.175
  
class Peak(object):
  def __init__(self, hs, ks, I, sigma=1):
    if len(hs) == len(ks):
      self.hs = np.array(hs, int)
      self.ks = np.array(ks, int)
      self.I = float(I)
      self.sigma = float(sigma)
    else:
      print("Dude, length of hs and ks are different")

 
def read_data_5_columns(filename="ripple_082-085.dat"):
    """Read a five-column ASCII file and parse each column into a python list.
    Lines starting with # will be ignored, i.e., # signals a comment line.
    
    filename: input file name
    
    The input file must be formatted as "h k q I", where 
    h, k: ripple main and side peak index
    q: magnitude of scattering vector, q
    I: observed intensity
    sigma: uncertainty in intensity
    For example, an input file should look like:
    
    # Example 1
    # =========
    # Comment goes here
    # Another comment line
    h  k      q      I  sigma
    1 -1  0.107   78.9    8.1
    1  0  0.100  100.0   10.0
    2  0  0.200   45.6    6.9
    2  1  0.205   56.7    8.0
    
    # Example 2 (include combined peaks)
    # In this case, separate indices by comma without any space.
    # The example shows a case in which different k's are combined.
    # ==========================================================
        h      k         q      I  sigma
      1,1   -1,0  combined  183.4      1
        1      1    0.1241   43.8      1
    2,2,2 -1,0,1  combined  180.0      1
    
    # Example 3 (combine different h's)
    # =======================================
          h         k         q      I  sigma
        1,1      -1,0  combined  183.4      1
    1,2,2,2  1,-1,0,1  combined  223.8      1
    """
    # Process comment and header lines
    fileobj = open(filename, 'r')
    while True:
        s = fileobj.readline()
        if s.startswith('#'):
            print(s),
            continue
        elif s.startswith('h'):
            break
        else:
            print("Any comments (including an empty line) should start with #.")
            print("Please fix your input file.")
            sys.exit(1)
    print("")
    
    # Go through data points  
    hl = []; kl = []; ql = []; Il = []; sl =[]; combl = []
    lines = fileobj.readlines()
    counter = 1
    for line in lines:
        # This ignores an empty line
        line = line.rstrip()
        if not line: 
            continue
        h, k, q, I, s = line.split()
        h = map(int, h.split(','))
        k = map(int, k.split(','))
        I = float(I)
        s = float(s)
        if len(h) == 1 and len(k) == 1:
            q = float(q)
            hl.extend(h); kl.extend(k); ql.append(q); Il.append(I); sl.append(s)
        elif len(h) != len(k):
            print("Please check line {0} to make sure that h and k ".format(counter))
            print("columns have the same number of items. For example,")
            print("1,1,1 -1,0,1 to combine (1,-1), (1,0), and (1,1) peaks.")
            sys.exit(1)
        else:
            combl.append((tuple(h), tuple(k), I, s))
        counter += 1
      
    return hl, kl, ql, Il, sl, combl


def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

###############################################################################
class BaseRipple(object):
    """The base ripple class, which will be inherited by contour subclasses, which
    in turn will be inherited by transbilayer subclasses. This base class mainly
    deals with the ripple lattice parameters, namely, D, lambda_r, and gamma.

    h, k: ripple phase main and side peak index
    qx, qz, q: scattering vector (qx is approx. equal to qr)
    I: observed intensity of each peak, before geometric correction
    D, lambda_r, gamma: D-spacing, ripple wavelength, and oblique angle
    sigma: uncertainties in intensity (equal to sqrt(I) = F)
    mask: mask out peaks that are set to False. Masked peak will be excluded from 
          the nonlinear fit. Only work for individual peaks, not combined ones.
    comb: a list of combined peaks. Each element in the list is a tuple,
          ((list of h index), (list of k index), F), where lists of index tells
          which peaks are combined to give the value of F. 
    """
    def __init__(self, h, k, q=None, I=None, sigma=None, D=58, lambda_r=140, 
                 gamma=1.7):
        self.h = np.array(h, int)
        self.k = np.array(k, int)
        self.q = np.array(q, float)
        self.I = np.array(I, float)
        self.F = sqrt(np.array(I, float))
        self.phase = np.ones(self.F.size)
        self.sigma = np.array(sigma, float)
        self.latt_par = Parameters()
        self.latt_par.add('D', value=D, vary=True)
        self.latt_par.add('lambda_r', value=lambda_r, vary=True)
        self.latt_par.add('gamma', value=gamma, vary=True)
        self.mask = np.ones(self.h.size, dtype=bool)
        self._set_qxqz()
        self.edm = ElectronDensityMap()
        self.edm.qx = self.qx
        self.edm.qz = self.qz
        self.edm.F = self.phase * self.F
    
    def show_mask(self):
        """Show the mask array."""
        print(" h  k  mask")
        for a, b, c in zip(self.h, self.k, self.mask):
            print("{0: 1d} {1: 1d}  {2:s}".format(a, b, c))
      
    def set_mask(self, h, k, value):
        """Set a mask element to True/False"""
        self.mask[(self.h==h)&(self.k==k)] = value
      
    def get_mask(self, h, k):
        """Get a mask element"""
        return self.mask[(self.h==h)&(self.k==k)]
    
    def _set_phase(self):
        """Model method needs to be defined in the subclasses"""
        self.phase = np.sign(self._model_F())
        self.phase = self.phase.astype(int)
        
    def flip_phase(self, h, k):
        """Flip the phase factor of the (h,k) order"""
        self.phase[(self.h==h)&(self.k==k)] *= -1
        
    def get_phase(self, h, k):
        """Get the phase factor of the (h,k) order"""
        return self.phase[(self.h==h)&(self.k==k)]
       
    def apply_Lorentz_correction(self, I):
        """Apply the Lorentz correction to the input intensity and return it. 
        
        I: observed intensity, which will be Lorentz corrected.
        """
        global wavelength
        ret = np.array(I)
        ret[self.k==0] = ret[self.k==0] * self.qz[self.k==0]
        ret[self.k!=0] = ret[self.k!=0] / wavelength / self.qz[self.k!=0] * 4 * \
                       np.pi * np.pi * np.absolute(self.qx[self.k!=0])
        return ret
    
    def apply_Lorentz_factor(self, I):
        """Apply the Lorentz factor to the input intensity and return it.
        I: form factor squared."""
        global wavelength
        ret = np.array(I)
        ret[self.k==0] = ret[self.k==0] / self.qz[self.k==0]
        ret[self.k!=0] = ret[self.k!=0] * wavelength * self.qz[self.k!=0] / 4 / \
                     np.pi / np.pi / np.absolute(self.qx[self.k!=0])
        return ret
  
    def report_model_F(self):
        """Show the model form factor along with the experimental one, which is 
        normalized at (h=1,k=0).
        """
        model_F = self._model_F()
        exp_F = self.F
        # exp_F is normalized at (h=1,k=0) peak
        expF_10 = exp_F[(self.h==1)&(self.k==0)]
        exp_F = exp_F / expF_10 * 100
        model_F = model_F / expF_10 * 100
        print(" h  k      q   F_exp F_model")
        for a, b, c, d, e in zip(self.h, self.k, self.q, exp_F, model_F):
            print("{0: 1d} {1: 1d} {2: .3f} {3: 7.2f} {4: 7.2f}".format(a, b, c, d, e))
      
    def report_model_I(self):
        """Show the model observed intensity along with the experimental I. Need a model
        method to call, which should be implemented in a derived class.
        """
        chi_square = ((self._model_intrinsic_I()-self.I) / self.sigma) ** 2
        print(" h  k     qx     qz      q      model          I     sigma     chi^2")
        for a, b, c, d, e, f, g, h, i in zip(self.h, self.k, self.qx, self.qz, self.q, self._model_intrinsic_I(), self.I, self.sigma, chi_square):
            print("{0: 1d} {1: 1d} {2: .3f} {3: .3f} {4: .3f} {5: 10.0f} {6: 10.0f} {7: 9.0f} {8: 9.0f}".format(a, b, c, d, e, f, g, h, i))
        print("\nTotal chi^2 = {0: .0f}".format(np.sum(chi_square)))
  
    def report_calc_lattice(self):
        """Show the calculated (fitted) q values for each peak along with 
        the input data."""
        print(" h  k  q_obs q_calc")
        q_calc = np.sqrt(self.q_square())
        for a, b, c, d in zip(self.h, self.k, self.q, q_calc):
            print("{0: 1d} {1: 1d} {2: .3f} {3: .3f}".format(a, b, c, d))
  
    def report_lattice(self):
        """Report the best fit lattice parameters"""
        lmfit.report_fit(self.latt_par)
        print("chisqr = {0:.3f}".format(self.lattice.chisqr))
          
    def fit_lattice(self):
        """Start a non-linear least squared fit for lattice parameters."""
        self.lattice = minimize(self._residual_lattice, self.latt_par)    

    def _residual_lattice(self, params):
        """params is a dummy variable, necessary for lmfit.minimize to work."""
        model = np.sqrt(self.q_square())
        data = np.absolute(self.q)
        return (model[self.mask] -data[self.mask])
       
    def q_square(self):
        """Return q^2 = qx^2 + qz^2 using the values of lambda_r, D, and gamma."""
        return self.q_x()**2 + self.q_z()**2    
  
    def q_x(self):
        """Return qx value in the ripple phase LAXS."""
        lambda_r = self.latt_par['lambda_r'].value 
        return 2*np.pi*self.k/lambda_r
  
    def q_z(self):
        """Return qz value in the ripple phase LAXS."""
        D = self.latt_par['D'].value
        lambda_r = self.latt_par['lambda_r'].value
        gamma = self.latt_par['gamma'].value
        return 2*np.pi*(self.h/D - self.k/lambda_r/np.tan(gamma))
  
    def _set_qxqz(self, h=None, k=None):
        """Set qx and qz arrays with the value of D, lambda_r, and gamma, 
        given lists of h and k indicies.
        """
        D = self.latt_par['D'].value
        lambda_r = self.latt_par['lambda_r'].value
        gamma = self.latt_par['gamma'].value
        if h is None:
          h = self.h
        if k is None:
          k = self.k
        self.qx = 2*np.pi*k/lambda_r
        self.qz = 2*np.pi*(h/D - k/lambda_r/np.tan(gamma))
    
    def report_edp(self):
        """Report best fit parameters related to electron density profile."""
        lmfit.report_fit(self.edp_par)
        print("chisqr = {0:.3f}".format(self.edp.chisqr))

    def fit_edp(self):
        """Start a non-linear least squared fit for electron density profile."""
        self._set_qxqz()
        self.edp = minimize(self._residual_edp, self.edp_par)
        self._set_phase()
    
    def _residual_edp(self, params):
        """Return the individual residuals."""
        model = self._model_intrinsic_I()
        data = self.I
        sigma = self.sigma
        return (data-model) / sigma 
          
    def _model_observed_I(self):
        """Apply the Lorentz factors to |model F|^2. Then return the properly 
        scaled, calculated observed intensity. After a fit is performed, this
        function returns the best fit.
        """
        global wavelength
        common_scale = self.edp_par['common_scale'].value
        I = self._model_F()**2
        I[self.k==0] = I[self.k==0] / self.qz[self.k==0]
        I[self.k!=0] = I[self.k!=0] * wavelength * self.qz[self.k!=0] / 4 / \
                       np.pi / np.pi / np.absolute(self.qx[self.k!=0])
        #I_10 = I[(self.h==1)&(self.k==0)]
        #I = I / I_10 * 10000 * common_scale
        I = I * common_scale
        return I
    
    def _model_intrinsic_I(self):
        """Intrinsic intensity, which is simply |F|^2"""
        I = self._model_F()**2
        return I
    
    def _model_F(self):
        """Return the model form factor. The return object is a numpy array with
        its length equal to the length of qx.
        """ 
        model = self.F_model()
        return model
  
    def export_model_F(self, outfilename="best_fit_F.dat"):
        """Export the best fit model form factor as an ASCII file consisting of 
        four columns, h, k, F_exp, F_model.
        
        outfilename: output file name
        """
        model_F = self._model_F()
        exp_F = self.F
        # exp_F is normalized at (h=1,k=0) peak
        expF_10 = exp_F[(self.h==1)&(self.k==0)]
        exp_F = exp_F / expF_10 * 100
        model_F = model_F / expF_10 * 100      
        with open(outfilename, 'w') as f:
            f.write(" h  k      q   F_exp F_model\n")
            for a, b, c, d, e in zip(self.h, self.k, self.q, exp_F, model_F):     
                f.write("{0: 1d} {1: 1d} {2: .3f} {3: 7.2f} {4: 7.2f}\n".format(a, b, c, d, e))

    def export_model_I(self, outfilename="best_fit_I.dat"):
        """Export the observed intensity calculated from a model as an ASCII file 
        consisting of nine columns.
        """
        chi_square = ((self._model_intrinsic_I()-self.I) / self.sigma) ** 2
        with open(outfilename, 'w') as ff:
            ff.write(" h  k     qx     qz      q    I_model   I_exp  sigma   chi^2\n")
            for a, b, c, d, e, f, g, h, i in zip(self.h, self.k, self.qx, 
                                                 self.qz, self.q, self._model_observed_I(), 
                                                 self.I, self.sigma, chi_square):
                ff.write("{0: 1d} {1: 1d} {2: .3f} {3: .3f} {4: .3f} {5: 8.0f} {6: 8.0f} \
                          {7: 5.0f} {8: 9.0f}\n".format(a, b, c, d, e, f, g, h, i))
            ff.write("\nTotal chi^2 = {0: .0f}".format(np.sum(chi_square)))

    def export_params(self, outfilename="params.txt"):
        with open(outfilename, 'w') as f:
            f.write("chisqr={0: f}\n".format(self.edp.chisqr))
            f.write("D={0: f}\n".format(self.latt_par['D'].value))
            f.write("lambda_r={0: f}\n".format(self.latt_par['lambda_r'].value))
            f.write("gamma={0: f}\n\n".format(self.latt_par['gamma'].value))      
            f.write(lmfit.fit_report(self.edp_par))
 
    def export_phases(self, filename):
        """
        Export an ASCII file containing the phases that reflect whatever phases
        the object has, which are not necessarily the same as
        the phases predicted by the model.
        """
        with open(filename, 'w') as f:
            f.write("h k phase\n")
            for a, b, c in zip(self.h, self.k, self.phase):
                f.write("{0: 1d} {1: 1d} {2: 1d}\n".format(a, b, c))
      
    def export_EDM(self, filename="EDM.dat", xmin=-150, xmax=150, zmin=-100,
                   zmax=100, N=201):
        """Export the Fourier-reconstructed 2D electron density map (EDM) as 
        an ASCII file consisting of three columns, x, z, and ED. 
        Calculate ED at N points along x and N points along z. The units are 
        in Angstrom.
        """
        self.edm.plot_EDM(xmin, xmax, zmin, zmax, N, self.phase*self.F, filename)
      
    def plot_EDM(self, xmin=-150, xmax=150, zmin=-100, zmax=100, N=201):
        """Plot an experimental 2D electron density map. Calculate
        EDM on an N by N grid. The units are in Angstrom.
        """    
        self.edm.plot_EDM(xmin, xmax, zmin, zmax, N, self.phase*self.F)

    def plot_model_EDM(self, xmin=-150, xmax=150, zmin=-100, zmax=100, N=201):
        """Plot a model 2D electron density map. Calculate
        EDM on an N by N grid. The units are in Angstrom.
        """    
        self.edm.plot_EDM(xmin, xmax, zmin, zmax, N, self._model_F())

    def export_EDP(self, filename="EDP.dat", center=(0,0), angle=-10, 
                   length=60, stepsize=0.5):
        self.edm.plot_EDP_angle(center, angle, length, stepsize, self.phase*self.F, 
                                filename)
  
    def export_model_EDP(self, filename="model_EDP.dat", center=(0.0), angle=-10, 
                         length=60, stepsize=0.5):
        self.edm.plot_EDP_angle(center, angle, length, stepsize, self._model_F(), 
                                filename)
  
    def plot_EDP(self, center=(0,0), angle=-10, length=60, stepsize=0.5):
        """Plot EDP along a line making an angle in degrees with respect to 
        the stacking z direction. Positive angle means the line is tilted in 
        the CW direction from the z-axis in the x-z plane. center specifies 
        about what position the plot is made. 
        Call plt.show() or plt.show(block=False) to actually display the plot.
	    
        Parameters
        ==========
        length: the total length in Angstrom, symmetric about center.
        stepsize: the step size in Angstrom in x.
          
        Example 
        =======
        plot_EDP(center=(0,0), angle=-10, length=60, stepsize=1)
          
        will plot the EDP about (x,z)=(0,0), along a line making 10 degrees
        in CCW from the z-axis, with +/-30 Angstrom above and below the
        center, calculated every Angstrom.
        """
        return self.edm.plot_EDP_angle(center, angle, length, stepsize, self.phase*self.F)
    
    def plot_EDP_between_two_points(self, start, end, N):
        """Plot an experimental EDP along a line connecting two points 
        secified by start and end, on N points. 
        """  
        return self.edm.plot_EDP_endpoints(start, end, N, self.phase*self.F)

    def export_EDP_between_two_points(self, filename, start, end, N):
        """Export an experimental EDP along a line connecting two points 
        secified by start and end, on N points, as an ASCII file.
        """
        self.edm.plot_EDP_endpoints(start, end, N, self.phase*self.F, filename)
    
    def get_EDP_between_two_points(self, start, end, N):
        self.edm.F = self.phase * self.F
        return self.edm.get_EDP_endpoints(start, end, N)
        
    def plot_model_EDP(self, center=(0,0), angle=0, length=60, stepsize=0.5):
        self.edm.plot_EDP_angle(center, angle, length, stepsize, self._model_F())
        
    def export_headgroup_positions(self, filename="temp.txt"): 
        """Export an ASCII file containing headgroup positions in both lower
        and upper leaflets. The first column is for x, the second for
        lower leaflet, and the third for upper leaflet. This method assumes
        that headgroups have the maximum electron density.    
        """
        self.edm.F = self.phase * self.F
        self.edm.export_headgroup_positions(self.latt_par['lambda_r'].value, 
                                            self.latt_par['D'].value, 
                                            self.edp_par['A'].value, 
                                            self.edp_par['xM'].value, 
                                            filename)

    def export_methyl_positions(self, filename="temp.txt"):
        """Export an ASCII file containing teminal methyl positions. 
        This method assumes that terminal methyls have the minimum electron density.       
        """
        self.edm.F = self.phase * self.F
        self.edm.export_methyl_positions(self.latt_par['lambda_r'].value, 
                                         self.latt_par['D'].value, 
                                         self.edp_par['A'].value, 
                                         self.edp_par['xM'].value, 
                                         filename)
                                
                                
###############################################################################
class ElectronDensityMap(object):
    """Implements electron density map (EDM) related methods. Normally,
    this class will be used inside the BaseRipple class, and deals with
    anything related to an electron density map. 
    """      
    def __init__(self):
        pass
                   
    def plot_EDM(self, xmin, xmax, zmin, zmax, N, F, filename=None):
        X, Y, Z = self._calc_EDM(xmin, xmax, zmin, zmax, N, F)
        if filename is None:
            X.shape = (N, N)
            Y.shape = (N, N)
            Z.shape = (N, N)
            plt.figure()
            rotate_Z = ndimage.rotate(Z, 90)
            imgplot = plt.imshow(rotate_Z, extent=[-150,150,-100,100], cmap='gray')
            return imgplot            
        else:
            with open(filename, 'w') as f:
                f.write("x z ED\n")
                for x, y, z in zip(X, Y, Z):
                    f.write("{0: 3.1f} {1: 3.1f} {2: }\n".format(x, y, z))
                
    def _calc_EDM(self, xmin, xmax, zmin, zmax, N, F):
        """Fourier-reconstruct a 2D map of the electron density profile and return
        as Z on (X,Y) grids. Calculate EDP at N points along x and N points along z. 
        The units are in Angstrom.
    
        output: X, Y, Z, each being numpy array
        """
        rho_xz = []
        xgrid = np.linspace(xmin, xmax, num=N)
        zgrid = np.linspace(zmin, zmax, num=N)  
        for x in xgrid:
            for z in zgrid:
                tmp = F * np.cos(self.qx*x+self.qz*z)
                rho_xz.append([x, z, tmp.sum(axis=0)])
        rho_xz = np.array(rho_xz, float)  
        X, Y, Z= rho_xz[:,0], rho_xz[:,1], rho_xz[:,2]
        return X, Y, Z
    
    def get_EDP_endpoints(self, start, end, N):
        x0, z0 = start
        x1, z1 = end
        xpoints = np.linspace(x0, x1, N)
        zpoints = np.linspace(z0, z1, N)
        X, Z, DIST, EDP = self._calc_EDP(xpoints, zpoints, start)
        return X, Z, DIST, EDP
        
    
    def plot_EDP_endpoints(self, start, end, N, F, filename=None): 
        """Plot an experimental EDP along a line connecting start 
        and end, on N points. If filename is specified, export an
        ASCII file instead.
        """
        x0, z0 = start
        x1, z1 = end
        xpoints = np.linspace(x0, x1, N)
        zpoints = np.linspace(z0, z1, N)
        return self._plot_EDP(xpoints, zpoints, start, F, filename)
            
    def plot_EDP_angle(self, center, angle, length, stepsize, F, filename=None):
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
        return self._plot_EDP(xpoints, zpoints, center, F, filename)     
    
    def _plot_EDP(self, xarray, zarray, center, F, filename):
        self.F = F
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
        
    def _calc_EDP(self, xpoints, zpoints, (xM,z0)):
        rho = []
        for x, z in zip(xpoints, zpoints):
            tmp = self.F * np.cos(self.qx*x+self.qz*z)
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
class Sawtooth(BaseRipple):
  def __init__(self, h, k, q, I, sigma, D=57.8, lambda_r=145, gamma=1.71, 
               xM=100, A=20):
    super(Sawtooth, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma)
    self.edp_par = Parameters()
    self.edp_par.add('xM', value=xM, vary=True)
    self.edp_par.add('A', value=A, vary=True)
    self.edp_par.add('f1', value=1, vary=False)
    self.edp_par.add('f2', value=0, vary=False)
  
  def F_model(self):
    """Ripple form factor"""
    f1 = self.edp_par['f1'].value
    f2 = self.edp_par['f2'].value
    common_scale = self.edp_par['common_scale'].value
    #sec = (-1)**self.k * (lr-xM) * sin(self.k*pi-w)/(self.k*pi-w)/lr
    return common_scale * (self.FCmajor()*self.FTmajor() + 
            f1*self.FCminor()*self.FTminor() + 
            f2*self.FCkink()*self.FTkink()) 
    
  def FCmajor(self):
    xM = self.edp_par['xM'].value
    A = self.edp_par['A'].value
    lr = self.latt_par['lambda_r'].value
    w = 0.5 * (self.qx*xM + self.qz*A)
    return xM * np.sin(w) / lr / w

  def FCminor(self):
    xM = self.edp_par['xM'].value
    A = self.edp_par['A'].value
    lr = self.latt_par['lambda_r'].value
    w = 0.5 * (self.qx*xM + self.qz*A)
    arg1 = 0.5*self.qx*lr + w
    arg2 = 0.5*self.qx*lr - w    
    return (lr-xM) * np.cos(0.5*arg1) * np.sin(arg2) / lr / np.cos(0.5*arg2) / arg2 
    
  def FCkink(self):
    xM = self.edp_par['xM'].value
    A = self.edp_par['A'].value
    lr = self.latt_par['lambda_r'].value
    w = 0.5 * (self.qx*xM + self.qz*A)
    return 2*np.cos(w)/lr
  
  def where_in_sawtooth(self, x):
    xM = self.edp_par['xM'].value
    A = self.edp_par['A'].value
    lr = self.latt_par['lambda_r'].value    
    return where_in_sawtooth(np.array([x]), lr, A, xM)
    
def where_in_sawtooth(x, lambda_r, A, xM):
  # First, make sure x is numpy.array
  x = np.asarray(x)
  if x.ndim == 0:
    x.shape = 1
  # Bring x outside of unit cell to the inside 
  while (x<-lambda_r/2).any() or (x>lambda_r/2).any():
    x[x<-lambda_r/2] = x[x<-lambda_r/2] + lambda_r
    x[x>lambda_r/2] = x[x>lambda_r/2] - lambda_r
  # z contains sawtooth z value corresponding to input x      
  z = np.zeros(x.size)
  z[x<-xM/2] = -A * (x[x<-xM/2] + lambda_r/2) / (lambda_r - xM)
  z[x>xM/2] = -A * (x[x>xM/2] - lambda_r/2) / (lambda_r - xM)
  z[~((x<-xM/2)|(x>xM/2))] = A * x[~((x<-xM/2)|(x>xM/2))] / xM
  return z
    
###############################################################################
class SDF(Sawtooth):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20, 
               common_scale=20, R_HM=2, X_h=20, psi=0.087):
    super(SDF, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, xM, A)
    self.edp_par.add('common_scale', value=common_scale, vary=True, min=0)
    self.edp_par.add('R_HM', value=R_HM, vary=True, min=0)
    self.edp_par.add('X_h', value=X_h, vary=True, min=0)
    self.edp_par.add('psi', value=psi, vary=True)
  
  def FTmajor(self):
    R_HM = self.edp_par['R_HM'].value
    X_h = self.edp_par['X_h'].value
    psi = self.edp_par['psi'].value  
    arg = self.qz*X_h*np.cos(psi) - self.qx*X_h*np.sin(psi)
    return (R_HM*np.cos(arg) - 1)
  
  def FTminor(self):
    return self.FTmajor()
    
  def FTkink(self):
    return self.FTmajor()

###############################################################################
class MDF(SDF):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20, f1=1, f2=0, 
               common_scale=20, R_HM=2, X_h=20, psi=0.087):
    super(MDF, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, xM, A, 
                              common_scale, R_HM, X_h, psi)
    self.edp_par['f1'].value = f1
    self.edp_par['f2'].value = f2
    self.edp_par['f1'].vary = True
    self.edp_par['f2'].vary = True


###############################################################################
class S2G(Sawtooth):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20.27, 
               rho_H1=2.21, Z_H1=20.00, sigma_H1=3.33,
               rho_H2=2.22, Z_H2=22.22, sigma_H2=3.33,
               rho_M=1, sigma_M=3, psi=0.087, common_scale=0.1):
    super(S2G, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, xM, A)
    self.edp_par.add('rho_H1', value=rho_H1, vary=True, min=0)
    self.edp_par.add('Z_H1', value=Z_H1, vary=True, min=0, max=30)
    self.edp_par.add('sigma_H1', value=sigma_H1, vary=True, min=0)
    self.edp_par.add('rho_H2', value=rho_H2, vary=True, min=0)
    self.edp_par.add('Z_H2', value=Z_H2, vary=True, min=0, max=30)
    self.edp_par.add('sigma_H2', value=sigma_H2, vary=True, min=0, max=10)
    self.edp_par.add('rho_M', value=rho_M, vary=True, min=0)    
    self.edp_par.add('sigma_M', value=sigma_M, vary=True, min=0, max=10)   
    self.edp_par.add('psi', value=psi, vary=True)
    self.edp_par.add('common_scale', value=common_scale, vary=True, min=0)    
  
  def FTmajor(self):
    rho_H1 = self.edp_par['rho_H1'].value
    Z_H1 = self.edp_par['Z_H1'].value
    sigma_H1 = self.edp_par['sigma_H1'].value
    rho_H2 = self.edp_par['rho_H2'].value
    Z_H2 = self.edp_par['Z_H2'].value
    sigma_H2 = self.edp_par['sigma_H2'].value
    rho_M = self.edp_par['rho_M'].value
    sigma_M = self.edp_par['sigma_M'].value
    psi = self.edp_par['psi'].value  
    
    return self.F2G(rho_H1, Z_H1, sigma_H1, rho_H2, Z_H2, sigma_H2,
                    rho_M, sigma_M, psi)

  def FTminor(self):
    return self.FTmajor()
    
  def FTkink(self):
    return self.FTmajor()
        
  def F2G(self, rho_H1, Z_H1, sigma_H1, rho_H2, Z_H2, sigma_H2,
              rho_M, sigma_M, psi):
    """
    Form factor calculated from the 2G hybrid model
    """  
    # Make sure Z_H2 > Z_H1. If Z_H2 < Z_H1, swap them
    if Z_H1 > Z_H2:
      Z_H1, Z_H2 = Z_H2, Z_H1
      sigma_H1, sigma_H2 = sigma_H2, sigma_H1
      rho_H1, rho_H2 = rho_H2, rho_H1
    
    # Calculate the intermediate variables
    alpha = self.qz*cos(psi) - self.qx*sin(psi)
    Z_CH2 = Z_H1 - sigma_H1
    Z_W = Z_H2 + sigma_H2
    DeltaZ_H = Z_W - Z_CH2
    
    # Calculate the Gaussian part   
    FG = -rho_M*sigma_M * exp(-0.5*(alpha*sigma_M)**2)
    FG += 2*rho_H1*sigma_H1 * cos(alpha*Z_H1) * exp(-0.5*(alpha*sigma_H1)**2)
    FG += 2*rho_H2*sigma_H2 * cos(alpha*Z_H2) * exp(-0.5*(alpha*sigma_H2)**2)
    FG *= np.sqrt(2*pi)
    
    # Calculate the strip part
    FS = -2 * sin(alpha*Z_CH2) / alpha
    
    # Calculate the bridging part
    FB = 1 / (alpha + pi/DeltaZ_H)
    FB += 1 / (alpha - pi/DeltaZ_H)
    FB *= sin(alpha*Z_W) + sin(alpha*Z_CH2)
    FB *= 0.5
    FB -= (sin(alpha*Z_W)-sin(alpha*Z_CH2)) / alpha
               
    return (FG + FS + FB)
    

###############################################################################
class M2G(S2G):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20, f1=1, f2=0, 
               rho_H1=2.21, Z_H1=20.24, sigma_H1=3.33,
               rho_H2=2.22, Z_H2=20.22, sigma_H2=3.33, 
               rho_M=1, sigma_M=3, psi=0.087, common_scale=0.1):
    super(M2G, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, xM, A, 
                              rho_H1, Z_H1, sigma_H1, rho_H2, Z_H2, sigma_H2, 
                              rho_M, sigma_M, psi, common_scale)
    self.edp_par['f1'].value = f1
    self.edp_par['f2'].value = f2
    self.edp_par['f1'].vary = True
    self.edp_par['f2'].vary = True


###############################################################################
class S1G(Sawtooth):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20.27, 
               rho_H_major=2.21, rho_H_minor=2.21,
               Z_H_major=20.00, Z_H_minor=20.00,
               sigma_H_major=3.33, sigma_H_minor=3.33,
               rho_M_major=1, rho_M_minor=1,
               sigma_M_major=3, sigma_M_minor=3,
               psi_major=0.087, psi_minor=0.087,
               common_scale=0.1):
    super(S1G, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, xM, A)
    self.edp_par.add('rho_H_major', value=rho_H_major, vary=True)
    self.edp_par.add('rho_H_minor', value=rho_H_minor, vary=False)
    self.edp_par.add('Z_H_major', value=Z_H_major, vary=True, min=0, max=30)
    self.edp_par.add('Z_H_minor', value=Z_H_minor, vary=False, min=0, max=30)
    self.edp_par.add('sigma_H_major', value=sigma_H_major, vary=True, min=0)
    self.edp_par.add('sigma_H_minor', value=sigma_H_minor, vary=False, min=0)
    self.edp_par.add('rho_M_major', value=rho_M_major, vary=True)    
    self.edp_par.add('rho_M_minor', value=rho_M_minor, vary=False) 
    self.edp_par.add('sigma_M_major', value=sigma_M_major, vary=True, min=0)  
    self.edp_par.add('sigma_M_minor', value=sigma_M_minor, vary=False, min=0)
    self.edp_par.add('psi_major', value=psi_major, vary=True)
    self.edp_par.add('psi_minor', value=psi_minor, vary=False)
    self.edp_par.add('common_scale', value=common_scale, vary=True) 
    self.link_rho_H = True
    self.link_Z_H = True
    self.link_sigma_H = True
    self.link_rho_M = True
    self.link_sigma_M = True
    self.link_psi = True  
  
  def unpack_major(self):
    rho_H_major = self.edp_par['rho_H_major'].value
    Z_H_major = self.edp_par['Z_H_major'].value
    sigma_H_major = self.edp_par['sigma_H_major'].value
    rho_M_major = self.edp_par['rho_M_major'].value
    sigma_M_major = self.edp_par['sigma_M_major'].value
    psi_major = self.edp_par['psi_major'].value  
    
    return rho_H_major, Z_H_major, sigma_H_major, rho_M_major, sigma_M_major, psi_major
    
  def unpack_minor(self):
    rho_H_minor = self.edp_par['rho_H_minor'].value
    Z_H_minor = self.edp_par['Z_H_minor'].value
    sigma_H_minor = self.edp_par['sigma_H_minor'].value
    rho_M_minor = self.edp_par['rho_M_minor'].value
    sigma_M_minor = self.edp_par['sigma_M_minor'].value
    psi_minor = self.edp_par['psi_minor'].value  
    
    return rho_H_minor, Z_H_minor, sigma_H_minor, rho_M_minor, sigma_M_minor, psi_minor
       
  def FTmajor(self):
    rho_H, Z_H, sigma_H, rho_M, sigma_M, psi = self.unpack_major() 
    return self.F1G(rho_H, Z_H, sigma_H, rho_M, sigma_M, psi)
    
  def FTminor(self):
    rho_H, Z_H, sigma_H, rho_M, sigma_M, psi = self.unpack_major()
    if self.link_rho_H is True:
      self.edp_par['rho_H_minor'].value = rho_H
    if self.link_Z_H is True:
      self.edp_par['Z_H_minor'].value = Z_H
    if self.link_sigma_H is True:
      self.edp_par['sigma_H_minor'].value = sigma_H
    if self.link_rho_M is True:
      self.edp_par['rho_M_minor'].value = rho_M
    if self.link_sigma_M is True:
      self.edp_par['sigma_M_minor'].value = sigma_M
    if self.link_psi is True:
      self.edp_par['psi_minor'].value = psi
    
    rho_H, Z_H, sigma_H, rho_M, sigma_M, psi = self.unpack_minor()
    return self.F1G(rho_H, Z_H, sigma_H, rho_M, sigma_M, psi)
    
  def FTkink(self):
    rho_H, Z_H, sigma_H, rho_M, sigma_M, psi = self.unpack_major()
    #rho_H, Z_H, sigma_H, rho_M, sigma_M, psi = self.unpack_minor() 
    return self.F1G(rho_H, Z_H, sigma_H, rho_M, sigma_M, psi)
  
  def F1G(self, rho_H1, Z_H1, sigma_H1, rho_M, sigma_M, psi):     
    # Calculate the intermediate variables
    alpha = self.qz*cos(psi) - self.qx*sin(psi)
    Z_CH2 = Z_H1 - sigma_H1
    Z_W = Z_H1 + sigma_H1
    DeltaZ_H = Z_W - Z_CH2
    
    # Calculate the Gaussian part   
    FG = -rho_M*sigma_M * exp(-0.5*(alpha*sigma_M)**2)
    FG += 2*rho_H1*sigma_H1 * cos(alpha*Z_H1) * exp(-0.5*(alpha*sigma_H1)**2)
    FG *= np.sqrt(2*pi)
    
    # Calculate the strip part
    FS = -2 * sin(alpha*Z_CH2) / alpha
    
    # Calculate the bridging part
    FB = 1 / (alpha + pi/DeltaZ_H)
    FB += 1 / (alpha - pi/DeltaZ_H)
    FB *= sin(alpha*Z_W) + sin(alpha*Z_CH2)
    FB *= 0.5
    FB -= (sin(alpha*Z_W)-sin(alpha*Z_CH2)) / alpha
               
    return (FG + FS + FB)


###############################################################################
class M1G(S1G):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20, f1=1, f2=0, 
               rho_H_major=2.21, rho_H_minor=2.21,
               Z_H_major=20.00, Z_H_minor=20.00,
               sigma_H_major=3.33, sigma_H_minor=3.33,
               rho_M_major=1, rho_M_minor=1,
               sigma_M_major=3, sigma_M_minor=3,
               psi_major=0.087, psi_minor=0.087,
               common_scale=0.1):
    super(M1G, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, xM, A, 
                              rho_H_major, rho_H_minor,
                              Z_H_major, Z_H_minor,
                              sigma_H_major, sigma_H_minor, 
                              rho_M_major, rho_M_minor,
                              sigma_M_major, sigma_M_minor,
                              psi_major, psi_minor,
                              common_scale)
    self.edp_par['f1'].value = f1
    self.edp_par['f2'].value = f2
    self.edp_par['f1'].vary = True
    self.edp_par['f2'].vary = True

  def export_params(self, outfilename="params.txt"):
    with open(outfilename, 'w') as f:
      f.write("chisqr={0: f}\n".format(self.edp.chisqr))
      f.write("D={0: f}\n".format(self.latt_par['D'].value))
      f.write("lambda_r={0: f}\n".format(self.latt_par['lambda_r'].value))
      f.write("gamma={0: f}\n\n".format(self.latt_par['gamma'].value))      
      f.write(lmfit.fit_report(self.edp_par))
              
  def export_params2(self, outfilename="params.txt"):
    with open(outfilename, 'w') as f:
      f.write("D={0: f}\n".format(self.latt_par['D'].value))
      f.write("lambda_r={0: f}\n".format(self.latt_par['lambda_r'].value))
      f.write("gamma={0: f}\n".format(self.latt_par['gamma'].value))
      f.write("xM={0: f}, {1: b}\n".format(self.edp_par['xM'].value, self.edp_par['xM'].vary))
      f.write("A={0: f}, {1: b}\n".format(self.edp_par['A'].value, self.edp_par['A'].vary))
      f.write("f1={0: f}, {1: b}\n".format(self.edp_par['f1'].value, self.edp_par['f1'].vary))
      f.write("f2={0: f}, {1: b}\n".format(self.edp_par['f2'].value, self.edp_par['f2'].vary))
      f.write("rho_H_major={0: f}, {1: b}\n".format(self.edp_par['rho_H_major'].value, self.edp_par['rho_H_major'].vary))
      f.write("rho_H_minor={0: f}, {1: b}\n".format(self.edp_par['rho_H_minor'].value, self.edp_par['rho_H_minor'].vary))
      f.write("Z_H_major={0: f}, {1: b}\n".format(self.edp_par['Z_H_major'].value, self.edp_par['Z_H_major'].vary))
      f.write("Z_H_minor={0: f}, {1: b}\n".format(self.edp_par['Z_H_minor'].value, self.edp_par['Z_H_minor'].vary))
      f.write("sigma_H_major={0: f}, {1: b}\n".format(self.edp_par['sigma_H_major'].value, self.edp_par['sigma_H_major'].vary))
      f.write("sigma_H_minor={0: f}, {1: b}\n".format(self.edp_par['sigma_H_minor'].value, self.edp_par['sigma_H_minor'].vary))
      f.write("rho_M_major={0: f}, {1: b}\n".format(self.edp_par['rho_M_major'].value, self.edp_par['rho_M_major'].vary))
      f.write("rho_M_minor={0: f}, {1: b}\n".format(self.edp_par['rho_M_minor'].value, self.edp_par['rho_M_minor'].vary))
      f.write("sigma_M_major={0: f}, {1: b}\n".format(self.edp_par['sigma_M_major'].value, self.edp_par['sigma_M_major'].vary))
      f.write("sigma_M_minor={0: f}, {1: b}\n".format(self.edp_par['sigma_M_minor'].value, self.edp_par['sigma_M_minor'].vary))
      f.write("psi_major={0: f}, {1: b}\n".format(self.edp_par['psi_major'].value, self.edp_par['psi_major'].vary))
      f.write("psi_minor={0: f}, {1: b}\n".format(self.edp_par['psi_minor'].value, self.edp_par['psi_minor'].vary))
      f.write("common_scale={0: f}, {1: b}\n".format(self.edp_par['common_scale'].value, self.edp_par['common_scale'].vary))
   
   
###############################################################################
def flatness(x):
    """Compute the flatness of a series of points"""
    chisqr = ((x-x.mean())**2).sum() / (x.size-1)
    return chisqr
    
def most_flat_profile(rip):
    mydict = {}
    best = 10**10
    # Grab indices for h = 6 orders
    #index = np.where((rip.h==6)|(rip.h==7)|(rip.h==9))
    index = np.where((rip.h==5))
    # There are N possible combinations for the phase factors
    N = 2**len(index[0])
    # Basic idea: convert i to binary string s. This loop will go through
    # all N possible combinations. 
    for i in xrange(N):
        s = generate_binary_string(i, N-1)
        a = np.array(binary_string_to_list(s), int)
        # Change binary 0 to phase factor -1
        a[a==0] = -1
        # Replace a subset of the phase factors with a generated combination
        rip.phase[index] = a       
        X, Z, DIST, ED = rip.get_EDP_between_two_points(start=(-40,21), end=(40,38), N=161)
        if flatness(ED) < best:
            best = flatness(ED)
            best_array = a
        mydict[s] = flatness(ED)
        #mylist.append(flatness(ED))
    rip.phase[index] = best_array
    return mydict, best_array, best
    
def binary_string_to_list(s):
    mylist = []
    for i in s:
        mylist.append(int(i))
    return mylist
    
def generate_binary_string(n, N):
    ret = "{0:b}".format(n)
    while len(ret) < len("{0:b}".format(N)):
        ret = '0' + ret
    return ret
  
def symmetry(x):
    center = np.argmin(x)
    for i in xrange(1,301):
        left = x[center-i]
        right = x[center+i]
    
    
def most_symmetric_profile(rip):
    pass
