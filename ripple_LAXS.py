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
  

class Base_Ripple:
  """ 
  The base ripple class 
  """
  def __init__(self, h, k, q, F):
    """
    
    """
    self.h = np.array(h, int)
    self.k = np.array(k, int)
    self.q = np.array(q, float)
    self.F = np.array(F, float)
    self.params = Parameters()
    self.params.add('D', value=58, vary=True)
    self.params.add('lambda_r', value=140, vary=True)
    self.params.add('gamma', value=1.7, vary=True)
        
  def fit_lattice(self):
    x = np.array([self.h, self.k])
    data = self.q * self.q
    self.lattice = minimize(self.residual_lattice, self.params)    

  def residual_lattice(self, params):
    D = params['D'].value
    lambda_r = params['lambda_r'].value
    gamma = params['gamma'].value
    
    model = self.model(D, lambda_r, gamma)
    data = self.q**2
    return (model -data)
  
  def update(self):
    self.qx = q_x(self.k, self.lambda_r)
    self.qz = q_z(self.h, self.k, self.D, self.lambda_r, self.gamma)
    
  def model(self, D, lambda_r, gamma):
 
  def q_square(h, k, D, l, g):
  """ 
  Return q = qx^2 + qz^2 in the ripple phase LAXS.
  
  h: h index
  k: k index
  D: D-spacing 
  l: lambda_r, the ripple wavelength
  g: gamma angle of the unit cell
  """
  h = self.h
  return q_x(k, l)*q_x(k, l) + q_z(h, k, D, l, g)*q_z(h, k, D, l, g)    
  
  def q_x(k, lambda_r):
  """
  Return qx value in the ripple phase LAXS.
  
  k: k index
  lambda_r: lambda_r, the ripple wavelength
  """
  return 2*np.pi*k/lambda_r
  
  def q_z(h, k, D, lambda_r, gamma):
  """
  Return qz value in the ripple phase LAXS.
  
  h: h index
  k: k index
  D: D-spacing 
  lambda_r: lambda_r, the ripple wavelength
  gamma: gamma angle of the unit cell
  """
  return 2*np.pi*(h/D - k/lambda_r/np.tan(gamma))
  
  
        
        












def omega(qx, qz, x0, A):
  """
  Return the intermediate variable, omega
  """
  return 0.5 * (qx*x0 + qz*A)
  

def SDF_model(params, qx, qz):
  """
  Return the simple delta function. The return object is a numpy array with
  its length equal to the length of qx.
  
  params: lmfit Parameter object, requires x0, A, rho_M, R_HM, X_h, psi, and
          lambda_r
  qx, qz: numpy array. These must be in the same length
  """
  x0 = params['x0'].value
  A = params['A'].value
  rho_M = params['rho_M'].value
  R_HM = params['R_HM'].value
  X_h = params['X_h'].value
  psi = params['psi'].value
  lambda_r = params['lambda_r'].value

  model = SDF_F_trans(qx, qz, X_h, psi, rho_M, R_HM) * \
          SDF_F_cont(qx, qz, x0, A, lambda_r)

  return model
  

def SDF_F_trans(qx, qz, X_h, psi, rho_M, R_HM):
  arg = qz*X_h*np.cos(psi) - qx*X_h*np.sin(psi)
  return rho_M * (R_HM*np.cos(arg) - 1)
  

def SDF_F_cont(qx, qz, x0, A, lambda_r):
  w = omega(qx, qz, x0, A)
  lr = lambda_r
  arg1 = 0.5*qx*lr + w
  arg2 = 0.5*qx*lr - w
  fir = x0 * np.sin(w) / lambda_r / w
  sec = (lr-x0) * np.cos(0.5*arg1) * np.sin(arg2) / lr / np.cos(0.5*arg2) / arg2 
  return (fir + sec)


def MDF_model(params, qx, qz):
  """
  
  """
  x0 = params['x0'].value
  A = params['A'].value
  rho_M = params['rho_M'].value
  R_HM = params['R_HM'].value
  X_h = params['X_h'].value
  psi = params['psi'].value
  lambda_r = params['lambda_r'].value
  f1 = params['f1'].value
  f2 = params['f2'].value
  
  model = MDF_F_trans(qx, qz, X_h, psi, rho_M, R_HM) * \
          MDF_F_cont(qx, qz, x0, A, lambda_r, f1, f2)
  
  return model

def MDF_F_trans(qx, qz, X_h, psi, rho_M, R_HM):
  """
  
  """
  return SDF_F_trans(qx, qz, X_h, psi, rho_M, R_HM)


def MDF_F_cont(qx, qz, x0, A, lambda_r, f1, f2):
  """
  
  """
  w = omega(qx, qz, x0, A)
  lr = lambda_r
  arg1 = 0.5*qx*lr + w
  arg2 = 0.5*qx*lr - w
  fir = x0 * np.sin(w) / lambda_r / w
  sec = (lr-x0) * np.cos(0.5*arg1) * np.sin(arg2) / lr / np.cos(0.5*arg2) / arg2 
  return (fir + f1*sec + 2*f2*np.cos(w)) 

    
# define objective function: returns the array to be minimized
def residual(params, model_type, qx, qz, data=None, sigma=None):
  """
  Return the individual residuals. If data=None, return the model. If 
  simga=None, return (data - model). Otherwise, return (data-model) / sigma.
  model is the form factor squared, that is, intensity. model_type can take
  'SDF' or 'MDF'.
  """
  
  if model_type == 'SDF':
    model = np.absolute(SDF_model(params, qx, qz))
  elif model_type == 'MDF':
    model = np.absolute(MDF_model(params, qx, qz))
    
  if data is None:
    return model
  if sigma is None:
    return (data - model) 
  return (data-model) / sigma



  
  
def Fourier_decomp(N=201, xmin=-100, xmax=100, zmin=-100, zmax=100):
  rho_xz = []
  xgrid = np.linspace(xmin, xmax, num=N)
  zgrid = np.linspace(zmin, zmax, num=N)
  for x in xgrid:
    for z in zgrid:
      tmp = phase * F * np.cos(qx*x+qz*z)
      rho_xz.append([x, z, tmp.sum(axis=0)])
  rho_xz = np.array(rho_xz, float)  
  X = rho_xz[:,0]
  Y = rho_xz[:,1]
  Z = rho_xz[:,2]
  X.shape = (N, N)
  Y.shape = (N, N)
  Z.shape = (N, N)
  plt.contourf(X, Y, Z)


def get_phase(params, qx, qz, model_type):
  if model_type == 'SDF':
    tmp = SDF_model(params, qx, qz)
  elif model_type == 'MDF':
    tmp = MDF_model(params, qx, qz)
    
  return np.sign(tmp)


def show_model():
  model = SDF_model(params, qx, qz)
  print " h  k      qx     qz      q   model       F"
  for a, b, c, d, e, f, g in zip(h, k, qx, qz, q, model, F):
    print("{0: 1d} {1: 1d}  {2: .3f} {3: .3f} {4: .3f} {5: 7.2f} {6: 7.2f}"
          .format(a, b, c, d, e, f, g))


def show_profile(start, end):
  """
  Show Fourier decomposed EDP along a line connecting the start and end points
  """
  rho = []
  x0, z0 = start
  x1, z1 = end
  N = 100 # number of points on which ED gets calculated
  xpoints = np.linspace(x0, x1, N)
  zpoints = np.linspace(z0, z1, N)
  for x, z in zip(xpoints, zpoints):
    tmp = phase * F * np.cos(qx*x+qz*z)
    dist = np.sqrt((x-x0)**2 + (z-z0)**2)
    rho.append([dist, tmp.sum(axis=0)])
  rho = np.array(rho, float)
  X = rho[:,0]
  Y = rho[:,1]
  plt.plot(X, Y)
  
  
if __name__ == "__main__":
  # read data to be fitted
  filename = "WackWebb2.dat"
  skip = 1
  infile = open(filename, 'r')
  h, k, q, F = read_data(infile, skip)
  h = np.array(h, int)
  k = np.array(k, int)
  q = np.array(q, float)
  F = np.array(F, float)
  
  # obtain lambda_r, D, and gamma using the Ripple class
  obj = Ripple(h, k, q, F)
  obj.fit()
  D = obj.D
  lambda_r = obj.lambda_r
  gamma = obj.gamma
  qx = q_x(k, lambda_r)
  qz = q_z(h, k, D, lambda_r, gamma)
  
  # Work on SDF
  params = Parameters()
  params.add('x0', value=105, vary=True)
  params.add('A', value=20.3, vary=True)
  params.add('rho_M', value=53.9, vary=True)
  params.add('R_HM', value=2.21, vary=True)
  params.add('X_h', value=20.2, vary=True)
  params.add('psi', value=0.0868, vary=True)
  params.add('lambda_r', value=lambda_r, vary=False)
  
  x = np.array(q, float)
  data = np.array(F, float)
  result = minimize(residual, params, args=('SDF', qx, qz, data))
  lmfit.report_fit(params)
  print("chisqr = {0:.3f}".format(result.chisqr))
  phase = get_phase(params, qx, qz, 'SDF')
  
  # Work on MDF
  params = Parameters()
  params.add('x0', value=103, vary=True)
  params.add('A', value=20, vary=True)
  params.add('rho_M', value=60, vary=True)
  params.add('R_HM', value=2, vary=True)
  params.add('X_h', value=20.4, vary=True)
  params.add('psi', value=0.157, vary=True)
  params.add('lambda_r', value=lambda_r, vary=False)
  params.add('f1', value=0.7, vary=True)
  params.add('f2', value=0, vary=True)
  
  result = minimize(residual, params, args=('MDF', qx, qz, data))
  lmfit.report_fit(params)
  print("chisqr = {0:.3f}".format(result.chisqr))
  phase = get_phase(params, qx, qz, 'MDF')
  
  
  # Work on S1G
  
  # Work on MDF
  
  # Work on M1G
  
  #Fourier_decomp(qx, qz, F, phase=None, N=201, xmin=-100, xmax=100, zmin=-100, zmax=100):
  
  # Optimization using the Ripple class
  #p = np.array([103, 18.6, 3, 20, 20.1, 0.0873])
  #obj.fit_SDF(par=p)
