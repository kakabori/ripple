import numpy as np
import matplotlib as ppl
import scipy.optimize

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
    hk = []
    qF = []
    for i in range(len(h)):
      hk.append([h[i], k[i]])
      qF.append([q[i], F[i]])
      
    self.hk = np.array(hk, int)
    self.qF = np.array(qF, float)
              
  def fit(self):
    x = self.hk
    y = self.qF[:,0] * self.qF[:,0]    
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
   
  x: an n-by-2 array, where the first column holds h values and the 
     second column k values; n is the number of data points
  D: D-spacing
  lambda_r: ripple wavelength 
  gamma: gamma angle in the unit cell """
  return q_square(x[:,0], x[:,1], D, lambda_r, gamma)


if __name__ == "__main__":
  filename = "WackWebb.dat"
  skip = 1
  infile = open(filename, 'r')
  h, k, q, F = read_data(infile, skip)
  obj = Ripple(h, k, q, F)
  obj.fit()
