import numpy as np
import matplotlib as ppl
import scipy.optimize
import matplotlib.pyplot as plt

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
  

class Ripple:
  """ 
  The ripple class for fitting LAXS data and getting electron density 
  profile. 
  """
  
  def __init__(self, h, k, q, F):
    """
    
    """
    self.h = np.array(h, int)
    self.k = np.array(k, int)
    self.q = np.array(q, float)
    self.F = np.array(F, float)

              
  def fit(self):
    x = np.array([self.h, self.k])
    y = self.q * self.q    
    result = scipy.optimize.curve_fit(func, x, y)    
    self.D = result[0][0]
    self.lambda_r = result[0][1]
    self.gamma = result[0][2]
    self.D_err = result[1][0,0]
    self.lambda_r_err = result[1][1,1]
    self.gamma_err = result[1][2,2]
    print "D:", self.D, "+/-", self.D_err
    print "lambda_r:", self.lambda_r, "+/-", self.lambda_r_err
    print "gamma:", self.gamma, "+/-", self.gamma_err


  def update(self):
    self.qx = q_x(self.k, self.lambda_r)
    self.qz = q_z(self.h, self.k, self.D, self.lambda_r, self.gamma)


  def Fourier_decomp(self, N=201, xmin=-100, xmax=100, zmin=-100, zmax=100):
    self.update()
    self.phase = np.array([-1,  1, -1, -1, -1, 1, 1, -1, -1, -1, 
                            1, -1,  1, -1,  1, 1, 1,  1, -1,  1])
    rho_xz = []
    xgrid = np.linspace(xmin, xmax, num=N)
    zgrid = np.linspace(zmin, zmax, num=N)
    for x in xgrid:
      for z in zgrid:
        tmp = self.phase * self.F * np.cos(self.qx*x+self.qz*z)
        rho_xz.append([x, z, tmp.sum(axis=0)])
    self.rho_xz = np.array(rho_xz, float)  
    X = self.rho_xz[:,0]
    Y = self.rho_xz[:,1]
    Z = self.rho_xz[:,2]
    X.shape = (N, N)
    Y.shape = (N, N)
    Z.shape = (N, N)
    plt.contourf(X, Y, Z)
    
  
  def fit_SDF(self, par=None, error=None):
    """
    Run this method only after running the fit method.
    
    par: Initial guess for the parameters. If None, the initial values will
         all be 1.
    error: sigma, which will be used as relative weights.
    """
    self.update()
    x = np.array([self.k, self.qx, self.qz])
    y = self.F * self.F    
    result = scipy.optimize.curve_fit(lambda x, x0, A, rho_M, R_HM, X_h, psi: SDF(x, x0, A, self.lambda_r, rho_M, R_HM, X_h, psi)**2, x, y, p0=par, sigma=error)    
    
    print result
    
    self.x0 = result[0][0]
    self.A = result[0][1]
    self.rho_M = result[0][2]
    self.R_HM = result[0][3]
    self.X_h = result[0][4]
    self.psi = result[0][5]
    
    self.x0_err = result[1][0,0]
    self.A_err = result[1][1,1]
    self.rho_M_err = result[1][2,2]
    self.R_HM_err = result[1][3,3]
    self.X_h_err = result[1][4,4]
    self.psi_err = result[1][5,5]
    
    print "x0:", self.x0, "+/-", self.x0_err
    print "A:", self.A, "+/-", self.A_err
    print "rho_M:", self.rho_M, "+/-", self.rho_M_err
    print "R_HM:", self.R_HM, "+/-", self.R_HM_err
    print "X_h:", self.X_h, "+/-", self.X_h_err
    print "psi:", self.psi, "+/-", self.psi_err
    
  def SDF(self):
    """
    Run this after fit_SDF
    """
    ret = F_T(self.qx, self.qz, self.rho_M, self.R_HM, self.X_h, self.psi) * \
          F_C(self.k, self.qx, self.qz, self.x0, self.A, self.lambda_r)
    return ret
    

def SDF(x, x0, A, lambda_r, rho_M, R_HM, X_h, psi):
  """
  Return simple delta function model.
  
  x: n-by-3 array. x(:,0) holds index k, x(:,1) qx, and x(:,2) qz values.
  """
  return F_T(x[1], x[2], rho_M, R_HM, X_h, psi) * \
         F_C(x[0], x[1], x[2], x0, A, lambda_r)


def F_T(qx, qz, rho_M, R_HM, X_h, psi):
  """
  Return the transbilayer part of the form factor
  """
  return rho_M * (R_HM*np.cos(qz*X_h*np.cos(psi)-qx*X_h*np.sin(psi)) - 1)


def F_C(k, qx, qz, x0, A, lambda_r):
  """
  Return the contour part of the form factor
  """ 
  return np.sin(omega(qx, qz, x0, A)) * \
         (np.pi*k*x0 - lambda_r*omega(qx, qz, x0, A)) / \
         omega(qx, qz, x0, A) / lambda_r / \
         (np.pi*k-omega(qx, qz, x0, A))
  
  
def omega(qx, qz, x0, A):
  """
  Return the intermediate variable, omega
  """
  return 0.5 * (qx*x0 + qz*A)
        
        
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


def q_square(h, k, D, l, g):
  """ 
  Return q = qx^2 + qz^2 in the ripple phase LAXS.
  
  h: h index
  k: k index
  D: D-spacing 
  l: lambda_r, the ripple wavelength
  g: gamma angle of the unit cell
  """
  return q_x(k, l)*q_x(k, l) + q_z(h, k, D, l, g)*q_z(h, k, D, l, g)    


def func(x, D, lambda_r, gamma):
  """ 
  A wrapper for q_square function to be used in the scipy.optimize.curve_fit().
   
  x: a 2-by-n array, where the first row holds h values and the 
     second row k values; n is the number of data points
  D: D-spacing
  lambda_r: ripple wavelength 
  gamma: gamma angle in the unit cell """
  return q_square(x[0], x[1], D, lambda_r, gamma)


if __name__ == "__main__":
  filename = "WackWebb.dat"
  skip = 1
  infile = open(filename, 'r')
  h, k, q, F = read_data(infile, skip)
  obj = Ripple(h, k, q, F)
  obj.fit()
  p = np.array([103, 18.6, 3, 20, 20.1, 0.0873])
  obj.fit_SDF(par=p)
