import sys
import numpy as np
from numpy import pi, sin, cos, tan, exp
import matplotlib as ppl
import scipy.optimize
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit


class Peak(object):
  def __init__(self, hs, ks, I, sigma=1):
    if len(hs) == len(ks):
      self.hs = np.array(hs, int)
      self.ks = np.array(ks, int)
      self.I = float(I)
      self.sigma = float(sigma)
    else:
      print("Dude, length of hs and ks are different")
  
class CombinedPeaks(object):
  """
  input hs and ks should be a list, tuple, or array
  I and sigma are floating points
  """
  def __init__(self):
    self.peak_list = []
    
  def add_peak(self, hs, ks, I, sigma=1):
    p = Peak(hs, ks, I, sigma)
    self.peak_list.append(p)
    
  def get_all_hkIsigma(self):
    """
    Return hs, ks, I, sigma
    """
    ret1 = []
    ret2 = []
    ret3 = []
    ret4 = []
    for p in self.peak_list:
      np.append(ret1, p.hs)
      np.append(ret2, p.ks)
      np.append(ret3, p.I)
      np.append(ret4, p.sigma)
    return ret1, ret2, ret3, ret4
    
  def shrink(self, arr):
    """
    Shrink the input array, arr, according to the prescription given by
    peak_list
    """
    index = 0
    ret = []
    for p in self.peak_list:
      ret.append(np.sum(arr[index:index+len(p.hs)]))
      index = index + len(p.hs)
      
    return np.array(ret)

 
def read_data_5_columns(filename="ripple_082-085.dat"):
  """
  Read a five-column ASCII file and parse each column into a python list.
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
      print(s)
      continue
    elif s.startswith('h'):
      break
    else:
      print("Any comments (including an empty line) should start with #.")
      print("Please fix your input file.")
      sys.exit(1)
  
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
  """ 
  The base ripple class, which will be inherited by contour subclasses, which
  in turn will be inherited by transbilayer subclasses. This base class mainly
  deals with the ripple lattice parameters, namely, D, lambda_r, and gamma.
  Methods such as showing 1D and 2D edp are also implemented.
  
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
    self.sigma = np.array(sigma, float)
    self.latt_par = Parameters()
    self.latt_par.add('D', value=D, vary=True)
    self.latt_par.add('lambda_r', value=lambda_r, vary=True)
    self.latt_par.add('gamma', value=gamma, vary=True)
    self.mask = np.ones(self.h.size, dtype=bool)
    self.comb_peaks = CombinedPeaks()

  def set_combined_peaks(self, comb):
    """
    comb is a list of lists. Within each are four elements. The first elememt
    is a tuple of h indicies. The second is a tuple of k indicies. The third
    is the sum of intensity of the individual indexed peaks. The last is
    the uncertainty on the sum of intensity.
    
    format: [[(h1a, h1b, h1c), (k1a, k1b, k1c), I1, sigma1], 
             [(h2a, h2b), (k2a, k2b), I2, sigma2], ...
            ]
    , where (h1a,k1a), (h1b,k1b), etc. are the Miller indicies. The first list
    in the example indicates that three peaks were combined to give the sum
    of I1 with sigma1 error. The second list means two peaks were combined.
    """
    for c in comb:
      self.comb_peaks.add_peak(hs=c[0], ks=c[1], I=c[2], sigma=c[3])
  
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
    """
    model method needs to be defined in the subclasses
    """
    self.phase = np.sign(self._model_F())
    
  def plot_2D_edp(self, xmin=-100, xmax=100, zmin=-100, zmax=100, N=201):
    """
    Plot Fourier-reconstruct a 2D map of the electron density profile. Calculate
    EDP at N points along x and N points along z. The units are in Angstrom.
    """
    rho_xz = []
    xgrid = np.linspace(xmin, xmax, num=N)
    zgrid = np.linspace(zmin, zmax, num=N)
    for x in xgrid:
      for z in zgrid:
        tmp = self.phase * self.F * np.cos(self.qx*x+self.qz*z)
        rho_xz.append([x, z, tmp.sum(axis=0)])
    rho_xz = np.array(rho_xz, float)  
    X, Y, Z= rho_xz[:,0], rho_xz[:,1], rho_xz[:,2]
    #Y = rho_xz[:,1]
    #Z = rho_xz[:,2]
    X.shape = (N, N)
    Y.shape = (N, N)
    Z.shape = (N, N)
    plt.figure()
    plt.contourf(X, Y, Z)
    
  def plot_1D_edp(self, start=(-10,25), end=(30,-20), N=100):
    """
    Plot Fourier-reconstructed EDP along a line connecting the start and end points
    N: number of points on which ED gets calculated
    Call plt.show() to actually display the plot.
    """
    rho = []
    x0, z0 = start
    x1, z1 = end
    xpoints = np.linspace(x0, x1, N)
    zpoints = np.linspace(z0, z1, N)
    for x, z in zip(xpoints, zpoints):
      tmp = self.phase * self.F * np.cos(self.qx*x+self.qz*z)
      dist = np.sqrt((x-x0)**2 + (z-z0)**2)
      rho.append([dist, tmp.sum(axis=0)])
    rho = np.array(rho, float)
    X = rho[:,0]
    Y = rho[:,1]
    plt.figure()
    plt.plot(X, Y)
  
  def _calc_F_from_observed_I(self):
    """
    The name says it all. Only display individual peaks, not combined ones.
    """
    global wavelength
    ret = np.array(self.I)
    self._set_qxqz(self.h, self.k)
    ret[self.k==0] = ret[self.k==0] * self.qz[self.k==0]
    ret[self.k!=0] = ret[self.k!=0] / wavelength / self.qz[self.k!=0] * 4 * \
                   np.pi * np.pi * np.absolute(self.qx[self.k!=0])
    ret = np.sqrt(ret)
    F_10 = ret[(self.h==1)&(self.k==0)]
    ret = ret / F_10 * 100
    return ret
  
  def report_model_F(self):
    """
    Show the model form factor.
    """
    self._set_qxqz(self.h, self.k)
    F = self._calc_F_from_observed_I()
    print(" h  k      q   model       F")
    for a, b, c, d, e in zip(self.h, self.k, self.q, self._model_F(), F):
      print("{0: 1d} {1: 1d} {2: .3f} {3: 7.2f} {4: 7.2f}".format(a, b, c, d, e))
      
  def report_model_I(self):
    """
    Show the model observed intensity along with the experimental I. Need a model
    method to call, which should be implemented in a derived class.
    """
    self._set_qxqz(self.h, self.k)
    print(" h  k     qx     qz      q      model          I")
    for a, b, c, d, e, f, g in zip(self.h, self.k, self.qx, self.qz, self.q, 
                                   self._model_observed_I(), self.I):
      print("{0: 1d} {1: 1d} {2: .3f} {3: .3f} {4: .3f} {5: 10.3f} {6: 10.3f}"
            .format(a, b, c, d, e, f, g))
  
  def report_calc_lattice(self):
    """
    Show the calculated (fitted) q values for each peak along with the input
    data
    """
    print(" h  k  q_obs q_calc")
    q_calc = np.sqrt(self.calc_q_square())
    for a, b, c, d in zip(self.h, self.k, self.q, q_calc):
      print("{0: 1d} {1: 1d} {2: .3f} {3: .3f}".format(a, b, c, d))
  
  def report_lattice(self):
    """
    Report the best fit lattice parameters
    """
    lmfit.report_fit(self.latt_par)
    print("chisqr = {0:.3f}".format(self.lattice.chisqr))
          
  def fit_lattice(self):
    """
    Start a non-linear least squared fit for lattice parameters.
    """
    self.lattice = minimize(self._residual_lattice, self.latt_par)    

  def _residual_lattice(self, params):
    """
    params is a dummy variable, necessary for lmfit.minimize to work
    """
    model = np.sqrt(self.calc_q_square())
    data = np.absolute(self.q)
    return (model[self.mask] -data[self.mask])
       
  def calc_q_square(self):
    """ 
    Return q^2 = qx^2 + qz^2 using the values of lambda_r, D, and gamma.
    """
    return self._q_x()**2 + self._q_z()**2    
  
  def _q_x(self):
    """
    Return qx value in the ripple phase LAXS.

    """
    lambda_r = self.latt_par['lambda_r'].value 
    return 2*np.pi*self.k/lambda_r
  
  def _q_z(self):
    """
    Return qz value in the ripple phase LAXS.
    """
    D = self.latt_par['D'].value
    lambda_r = self.latt_par['lambda_r'].value
    gamma = self.latt_par['gamma'].value
    return 2*np.pi*(self.h/D - self.k/lambda_r/np.tan(gamma))
  
  def _set_qxqz(self, h, k):
    """
    Set qx and qz arrays with the value of D, lambda_r, and gamma, given
    lists of h and k indicies.
    """
    D = self.latt_par['D'].value
    lambda_r = self.latt_par['lambda_r'].value
    gamma = self.latt_par['gamma'].value
    self.qx = 2*np.pi*k/lambda_r
    self.qz = 2*np.pi*(h/D - k/lambda_r/np.tan(gamma))
    
  def report_edp(self):
    """
    Report best fit parameters related to electron density profile
    """
    lmfit.report_fit(self.edp_par)
    print("chisqr = {0:.3f}".format(self.edp.chisqr))
        
  def fit_edp(self):
    """
    Start a non-linear least squared fit for electron density profile
    """
    self.edp = minimize(self._residual_edp, self.edp_par)
    self._set_phase()
    
  def _residual_edp(self, params):
    """
    Return the individual residuals.  
    """
    h, k, I, s = self.comb_peaks.get_all_hkIsigma()
    self._set_qxqz(np.append(self.h, h), np.append(self.k, k))
    
    model = self._model_observed_I()
    
    #Split the model to indivial and combined peaks
    model_indiv = model[0:len(self.h)]
    model_comb = model[len(self.h):]
    
    # Apply the mask
    model_indiv = model_indiv[self.mask]
    
    # Combine peaks according to prescription given by CombinedPeaks object
    model_comb = self.comb_peaks.shrink(model_comb)
    
    model = np.append(model_indiv, model_comb)
    
    # Apply the mask to data points for individual peaks
    data = self.I[self.mask]
    sigma = self.sigma[self.mask]
    data = np.append(data, I)
    sigma = np.append(sigma, s)

    return (data-model) / sigma 
        
    # The following three lines do not reproduce Sun's results, which proves
    # that the fits were done through intensity, not form factor.
    #data = self.F
    #model = np.absolute(self._model())
    #return (data - model) 
  
  def _model_observed_I(self):
    """
    Return the geometrically corrected observed intensity predicted by the 
    model form factor
    """
    global wavelength
    I = self._model_F()**2
    I[self.k==0] = I[self.k==0] / self.qz[self.k==0]
    I[self.k!=0] = I[self.k!=0] * wavelength * self.qz[self.k!=0] / 4 / \
                   np.pi / np.pi / np.absolute(self.qx[self.k!=0])
    I_10 = I[(self.h==1)&(self.k==0)]
    I = I / I_10 * 10000
    return I
    
  def _model_F(self):
    """
    Return the model form factor. The return object is a numpy array with
    its length equal to the length of qx.
    """
    rho_M = self.edp_par['rho_M'].value
    model = self.F_trans() * self.F_cont()
    # get F(h=1,k=0), which is used for normalization
    # rho_M is a common scaling factor => F(h=1,k=0) = 100*rho_M
    F_10 = model[(self.h==1)&(self.k==0)]
    model = model / np.absolute(F_10) * 100 * rho_M
    return model
    
  def export_2D_edp(self, filename="2Dedp.dat", xmin=-100, xmax=100, 
                    zmin=-100, zmax=100, N=201):
    """
    Export the Fourier-reconstructed 2D electron density (ED) map as an ASCII 
    file consisting of three columns, x, z, and ED. 
    Calculate EDP at N points along x and N points along z. The units are in 
    Angstrom.
    """
    rho_xz = []
    xgrid = np.linspace(xmin, xmax, num=N)
    zgrid = np.linspace(zmin, zmax, num=N)
    for x in xgrid:
      for z in zgrid:
        tmp = self.phase * self.F * np.cos(self.qx*x+self.qz*z)
        rho_xz.append([x, z, tmp.sum(axis=0)])
    rho_xz = np.array(rho_xz, float)  
    X, Y, Z= rho_xz[:,0], rho_xz[:,1], rho_xz[:,2]
    with open(filename, 'w') as f:
      f.write("x z ED\n")
      for x, y, z in zip(X, Y, Z):
        f.write("{0: 3.1f} {1: 3.1f} {2: }\n".format(x, y, z))

  def export_1D_edp(self, filename="1Dedp.dat", start=(-10,25), end=(30,-20), 
                    N=100):
    """
    Export the Fourier-reconstructed EDP along a line connecting the start and 
    end points as an ASCII file consisting of four columns, x, z, distance from
    the start point, and EDP.
    
    filename: output file name
    start=(xmax,zmax)
    end=(xmax,zmax)
    N: number of points on which ED gets calculated
    """
    rho = []
    x0, z0 = start
    x1, z1 = end
    xpoints = np.linspace(x0, x1, N)
    zpoints = np.linspace(z0, z1, N)
    for x, z in zip(xpoints, zpoints):
      tmp = self.phase * self.F * np.cos(self.qx*x+self.qz*z)
      dist = np.sqrt((x-x0)**2 + (z-z0)**2)
      rho.append([x, z, dist, tmp.sum(axis=0)])
    rho = np.array(rho, float)
    X, Z, DIST, EDP = rho[:,0], rho[:,1], rho[:,2], rho[:,3]
    with open(filename, 'w') as f:
      f.write("x z dist ED\n")
      for x, z, dist, edp in zip(X, Z, DIST, EDP):
        f.write("{0: 3.1f} {1: 3.1f} {2: 3.1f} {3: }\n".format(x, z, dist, edp))
    

###############################################################################
class Sawtooth(BaseRipple):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               x0=100, A=20):
    super(Sawtooth, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma)
    self.edp_par = Parameters()
    self.edp_par.add('x0', value=x0, vary=True)
    self.edp_par.add('A', value=A, vary=True)
    self.edp_par.add('f1', value=1, vary=False)
    self.edp_par.add('f2', value=0, vary=False)
  
  def F_cont(self):
    """
    Contuour part of the ripple form factor
    """
    x0 = self.edp_par['x0'].value
    A = self.edp_par['A'].value
    f1 = self.edp_par['f1'].value
    f2 = self.edp_par['f2'].value
    lr = self.latt_par['lambda_r'].value
    w = 0.5 * (self.qx*x0 + self.qz*A)
    arg1 = 0.5*self.qx*lr + w
    arg2 = 0.5*self.qx*lr - w
    fir = x0 * np.sin(w) / lr / w
    sec = (lr-x0) * np.cos(0.5*arg1) * np.sin(arg2) / lr / np.cos(0.5*arg2) / arg2 
    #sec = (-1)**self.k * (lr-x0) * sin(self.k*pi-w)/(self.k*pi-w)/lr
    return (fir + f1*sec + 2*f2*np.cos(w)/lr) 


###############################################################################
class SDF(Sawtooth):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               x0=100, A=20, 
               rho_M=20, R_HM=2, X_h=20, psi=0.087):
    super(SDF, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, x0, A)
    self.edp_par.add('rho_M', value=rho_M, vary=True)
    self.edp_par.add('R_HM', value=R_HM, vary=True)
    self.edp_par.add('X_h', value=X_h, vary=True)
    self.edp_par.add('psi', value=psi, vary=True)
  
  def F_trans(self):
    """
    Transbilayer part of the ripple form factor
    """
    rho_M = self.edp_par['rho_M'].value
    R_HM = self.edp_par['R_HM'].value
    X_h = self.edp_par['X_h'].value
    psi = self.edp_par['psi'].value  
    arg = self.qz*X_h*np.cos(psi) - self.qx*X_h*np.sin(psi)
    return rho_M * (R_HM*np.cos(arg) - 1)
  

###############################################################################
class MDF(SDF):
  def __init__(self, h, k, q, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               x0=100, A=20, f1=1, f2=0, 
               rho_M=20, R_HM=2, X_h=20, psi=0.087):
    super(MDF, self).__init__(h, k, q, I, sigma, D, lambda_r, gamma, x0, A, 
                              rho_M, R_HM, X_h, psi)
    self.edp_par['f1'].value = f1
    self.edp_par['f2'].value = f2
    self.edp_par['f1'].vary = True
    self.edp_par['f2'].vary = True


###############################################################################
class SGF(Sawtooth):
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, x0=100, A=20.27, R_HM=2.21, X_H=20.24, sigma_H=3, 
               rho_M=20, sigma_M=3, psi=0.087):
    super(SGF, self).__init__(h, k, F, q, qx, qz, D, lambda_r, gamma, x0, A)
    self.edp_par.add('R_HM', value=R_HM, vary=True)
    self.edp_par.add('X_H', value=X_H, vary=True)
    self.edp_par.add('sigma_H', value=sigma_H, vary=True)
    self.edp_par.add('rho_M', value=rho_M, vary=True)
    self.edp_par.add('sigma_M', value=sigma_M, vary=True)
    self.edp_par.add('psi', value=psi, vary=True)
  
  def F_trans(self):
    """
    Transbilayer part of the ripple form factor
    """
    R_HM = self.edp_par['R_HM'].value
    X_H = self.edp_par['X_H'].value
    sigma_H = self.edp_par['sigma_H'].value
    rho_M = self.edp_par['rho_M'].value
    sigma_M = self.edp_par['sigma_M'].value
    psi = self.edp_par['psi'].value  
    th = np.cos(psi)*(self.qz - self.qx*np.tan(psi))
    first = R_HM*sigma_H*np.cos(X_H*th)*np.exp(-0.5*sigma_H**2*th**2)
    second = sigma_M*np.exp(-0.5*sigma_M**2*th**2)
    return np.sqrt(2*pi)*np.cos(psi)*rho_M*(first-second)  

      
###############################################################################
class MGF(SGF): 
  pass


###############################################################################
class S1G(Sawtooth):
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, x0=100, A=20.27, R_HM=2.21, X_H=20.24, sigma_H=3, 
               rho_M=20, sigma_M=3, psi=0.087, drho=0.1):
    super(S1G, self).__init__(h, k, F, q, qx, qz, D, lambda_r, gamma, x0, A)
    self.edp_par.add('R_HM', value=R_HM, vary=True)
    self.edp_par.add('X_H', value=X_H, vary=True)
    self.edp_par.add('sigma_H', value=sigma_H, vary=True)
    self.edp_par.add('rho_M', value=rho_M, vary=True)
    self.edp_par.add('sigma_M', value=sigma_M, vary=True)
    self.edp_par.add('psi', value=psi, vary=True)
    self.edp_par.add('delta_rho', value=drho, vary=True)
  
  def F_trans(self):
    """
    Transbilayer part of the ripple form factor
    """
    R_HM = self.edp_par['R_HM'].value
    Z_H = self.edp_par['X_H'].value
    sigma_H = self.edp_par['sigma_H'].value
    rho_M = self.edp_par['rho_M'].value
    sigma_M = self.edp_par['sigma_M'].value
    psi = self.edp_par['psi'].value  
    drho = self.edp_par['delta_rho'].value
    th = cos(psi)*(self.qz - self.qx*tan(psi))
    Fs = 2*drho*sin(self.qz*Z_H)*cos(self.qz*w/2)
    Fs = Fs * (-1/self.qz + self.qz*w*w/(self.qz**2*w*w-pi*pi))
    FG = sqrt(2*pi)*cos(psi*rho_M)
    FG = FG * (R_HM*sigma_H*cos(Z_H*th)*exp(-0.5*sigma_H**2*th*th)
               -sigma_M*exp(-0.5*sigma_M**2*th*th))
    return (Fs + FG)
    

###############################################################################
class M1G(S1G):
  pass
    










def F_C(h=1,k=0,D=57.94,lr=141.7,gamma=1.7174,x0=103,A=18.6):
  qx = 2*pi*k/lr
  qz = 2*pi*h/D - 2*pi*k/lr/tan(gamma)
  w = 0.5*(qx*x0+qz*A)
  arg1 = 0.5*qx*lr + w
  arg2 = 0.5*qx*lr - w
  return x0*sin(w)/w/lr+(lr-x0)/lr*cos(0.5*arg1)/cos(0.5*arg2)*sin(arg2)/arg2  

def F_C2(h=1,k=0,D=57.94,lr=141.7,gamma=1.7174,x0=103,A=18.6):
  qx = 2*pi*k/lr
  qz = 2*pi*h/D - 2*pi*k/lr/tan(gamma)
  w = 0.5*(qx*x0+qz*A)
  return sin(w)*(k*pi*x0/lr-w)/(w*(k*pi-w))

def F_T(h=1,k=0,D=57.94,lr=141.7,gamma=1.7174,rhom=51.38,rhm=2.2,xh=20.1,psi=5):
  psi = psi * pi / 180
  qx = 2*pi*k/lr
  qz = 2*pi*h/D - 2*pi*k/lr/tan(gamma)
  return rhom*(rhm*cos(qz*xh*cos(psi)-qx*xh*sin(psi)) - 1)





if __name__ == "__main__":
#  reproduce_WenjunSun_PNAS()

  # read data to be fitted
#  infilename = 'ripple_082-085.dat'
#  infilename = 'WackWebb2.dat'
#  infilename = 'ripple_082-085_mod1.dat'
  wavelength = 1.175
  infilename = 'intensity/ripple_082-085.dat'
  h, k, q, I, sigma, combined = read_data_5_columns(infilename)

  
  # Work on SDF
  sdf = SDF(h, k, q, I, sigma, D=57.54, lambda_r=142.88, gamma=1.706, 
            x0=102.8, A=21.1, rho_M=1.04, R_HM=2.52, X_h=19.7, psi=0.0823) 
  sdf.set_combined_peaks(combined)
#  sdf.set_mask(h=1, k=2, value=False)
#  sdf.set_mask(h=3, k=5, value=False)
#  sdf.set_mask(h=3, k=6, value=False)
#  sdf.set_mask(h=4, k=0, value=False)
#  sdf.fit_lattice()
#  sdf.fit_edp()
#  sdf.report_edp()
  
  # Work on MDF
#  mdf = MDF(h, k, F, q, qx=None, qz=None, D=58, lambda_r=141.7, gamma=1.7174, 
#            x0=103, A=18.6, f1=1, f2=0, rho_M=1, R_HM=2.2, 
#            X_h=20.1, psi=0.0872) 
#  mdf.fit_lattice()
#  mdf.fit_edp()
#  mdf.report_edp()  
