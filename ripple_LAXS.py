import numpy as np
import matplotlib as ppl
import scipy.optimize
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit

def read_data(fileobj, skip=0):
  """
  Read a four-column ASCII file and parse each column into a python list.
  
  fileobj: input file object
  skip: the number of header lines that will be skipped
  """
  # ignore the first skip lines
  for i in range(skip):
    fileobj.readline()
  h = []; k = []; q = []; F = []
  lines = fileobj.readlines()
  for line in lines:
    hval, kval, qval, Fval = line.split()
    h.append(hval); k.append(kval); q.append(qval); F.append(Fval) 
  return h, k, q, F
  

###############################################################################
class BaseRipple(object):
  """ 
  The base ripple class, which will be inherited by contour subclasses, which
  in turn will be inherited by transbilayer subclasses. This base class mainly
  deals with the ripple lattice parameters, namely, D, lambda_r, and gamma.
  Methods such as showing 1D and 2D edp are also implemented.
  """
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7):
    self.h = np.array(h, int)
    self.k = np.array(k, int)
    self.q = np.array(q, float)
    self.F = np.array(F, float)
    self.qx = np.array(qx, float)
    self.qz = np.array(qz, float)
    self.latt_par = Parameters()
    self.latt_par.add('D', value=D, vary=True)
    self.latt_par.add('lambda_r', value=lambda_r, vary=True)
    self.latt_par.add('gamma', value=gamma, vary=True)

  def _set_phase(self):
    """
    model method needs to be defined in the subclasses
    """
    self.phase = np.sign(self._model())
    
  def show_2D_edp(self, N=201, xmin=-100, xmax=100, zmin=-100, zmax=100):
    """
    Fourier-reconstruct a 2D map of the electron density profile. Calculate
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
    plt.contourf(X, Y, Z)
    
  def show_1D_edp(self, start=(-10,25), end=(30,-20), N=100):
    """
    Show Fourier-reconstructed EDP along a line connecting the start and end points
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
    plt.plot(X, Y)
    
  def report_model_F(self):
    """
    Show the model form factor along with the experimental |F|. Need a model
    method to call, which should be implemented in a derived class.
    """
    print(" h  k      qx     qz      q   model       F")
    for a, b, c, d, e, f, g in zip(self.h, self.k, self.qx, self.qz, self.q, 
                                   self._model(), self.F):
      print("{0: 1d} {1: 1d}  {2: .3f} {3: .3f} {4: .3f} {5: 7.2f} {6: 7.2f}"
            .format(a, b, c, d, e, f, g))
  
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
    self.qx = self._q_x()
    self.qz = self._q_z()

  def _residual_lattice(self, params):
    """
    params is a dummy variable, necessary for lmfit.minimize to work
    """
    model = self.model_lattice()
    data = self.q**2
    return (model -data)
    
  def model_lattice(self):
    D = self.latt_par['D'].value
    lambda_r = self.latt_par['lambda_r'].value
    gamma = self.latt_par['gamma'].value   
    return self._q_square()
       
  def _q_square(self):
    """ 
    Return q = qx^2 + qz^2 in the ripple phase LAXS.
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
    self.edp = minimize(self._residual, self.edp_par)
    self._set_phase()
    
  def _residual(self, params):
    """
    Return the individual residuals.  
    """
    data = self.F
    model = np.absolute(self._model())
    return (data - model) 
    #return (data-model) / sigma    
    
  def _model(self):
    """
    Return the model form factor. The return object is a numpy array with
    its length equal to the length of qx.
    """
    model = self.F_trans() * self.F_cont()
    return model
  

###############################################################################
class Sawtooth(BaseRipple):
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, x0=100, A=20):
    super(Sawtooth, self).__init__(h, k, F, q, qx, qz, D, lambda_r, gamma)
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
    w = self._omega()
    arg1 = 0.5*self.qx*lr + w
    arg2 = 0.5*self.qx*lr - w
    fir = x0 * np.sin(w) / lr / w
    sec = (lr-x0) * np.cos(0.5*arg1) * np.sin(arg2) / lr / np.cos(0.5*arg2) / arg2 
    return (fir + f1*sec + 2*f2*np.cos(w)) 
  
  def _omega(self):
    """
    Return the intermediate variable, omega
    """
    x0 = self.edp_par['x0'].value
    A = self.edp_par['A'].value
    return 0.5 * (self.qx*x0 + self.qz*A)


###############################################################################
class SDF(Sawtooth):
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, x0=100, A=20, rho_M=20, R_HM=2, X_h=20, 
               psi=0.087):
    super(SDF, self).__init__(h, k, F, q, qx, qz, D, lambda_r, gamma, x0, A)
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
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, x0=100, A=20, f1=1, f2=0, rho_M=20, R_HM=2, X_h=20, 
               psi=0.087):
    super(MDF, self).__init__(h, k, F, q, qx, qz, D, lambda_r, gamma, x0, A, 
                     rho_M, R_HM, X_h, psi)
    self.edp_par['f1'].value = f1
    self.edp_par['f2'].value = f2
    self.edp_par['f1'].vary = True
    self.edp_par['f2'].vary = True


###############################################################################
class SGF(Sawtooth):
  def __init__(self, h, k, F, q=None, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, x0=100, A=20, rho_M=20, R_HM=2, X_h=20, 
               psi=0.087):
    super(SDF, self).__init__(h, k, F, q, qx, qz, D, lambda_r, gamma, x0, A)
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
class MGF(SGF): 
  pass


###############################################################################
class S1G(Sawtooth):
  pass


###############################################################################
class M1G(S1G):
  pass
    










  
  
if __name__ == "__main__":
  # read data to be fitted
  infile = open("WackWebb2.dat", 'r')
  h, k, q, F = read_data(infile, skip=1)
  h = np.array(h, int)
  k = np.array(k, int)
  q = np.array(q, float)
  F = np.array(F, float) 
  sdf = SDF(h, k, F, q, qx=None, qz=None, D=58, lambda_r=140, gamma=1.7, x0=105, 
            A=20.3, rho_M=53.9, R_HM=2.21, X_h=20.2, psi=0.0868)   
  sdf.fit_lattice()
  #sdf.report_lattice()
  sdf.fit_edp()
  sdf.report_edp()
  
  # Work on MDF
  mdf = MDF(h, k, F, q, qx=None, qz=None, D=58, lambda_r=140, gamma=1.7, 
            x0=103, A=20, f1=0.7, f2=0, rho_M=60, R_HM=2, X_h=20.4, psi=0.157)
  mdf.fit_lattice()
  #mdf.report_lattice()
  mdf.fit_edp()
  mdf.report_edp()
  
  # Work on S1G
  
  # Work on MDF
  
  # Work on M1G
  
  #Fourier_decomp(qx, qz, F, phase=None, N=201, xmin=-100, xmax=100, zmin=-100, zmax=100):
  
  # Optimization using the Ripple class
  #p = np.array([103, 18.6, 3, 20, 20.1, 0.0873])
  #obj.fit_SDF(par=p)
