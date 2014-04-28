import numpy as np
from numpy import pi, sin, cos, tan, exp
import matplotlib as ppl
import scipy.optimize
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit

def read_data(fileobj, skip=0, data_type="F"):
  """
  Read a four-column ASCII file and parse each column into a python list.
  
  fileobj: input file object
  skip: the number of header lines that will be skipped
  data_type: whether input data are form factor or intensity. If intensity,
             the correction for oriented sample will be applied.
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
  def __init__(self, h, k, F, q, qx=None, qz=None, D=58, lambda_r=140, 
               gamma=1.7, I=None):
    self.h = np.array(h, int)
    self.k = np.array(k, int)
    self.q = np.array(q, float)
    self.F = np.array(F, float)
    self.latt_par = Parameters()
    self.latt_par.add('D', value=D, vary=True)
    self.latt_par.add('lambda_r', value=lambda_r, vary=True)
    self.latt_par.add('gamma', value=gamma, vary=True)
    if qx is None:
      self.qx = self._q_x()
    else:
      self.qx = np.array(qx, float)
    if qz is None:
      self.qz = self._q_z()
    else:
      self.qz = np.array(qz, float)
    if I is None:
      self.I = self.F**2
    else:
      self.I = I
    self.sigma = np.absolute(self.F)
      
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
    data = self.F**2
    model = np.absolute(self._model())**2
    return (data-model) / self.sigma  
    # The following three lines do not reproduce Sun's results, which proves
    # that the fits were done through intensity, not form factor.
    #data = self.F
    #model = np.absolute(self._model())
    #return (data - model) 
      
    
  def _model(self):
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
  

###############################################################################
class Sawtooth(BaseRipple):
  def __init__(self, h, k, F, q, qx=None, qz=None, D=58, lambda_r=140, 
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
    w = 0.5 * (self.qx*x0 + self.qz*A)
    arg1 = 0.5*self.qx*lr + w
    arg2 = 0.5*self.qx*lr - w
    fir = x0 * np.sin(w) / lr / w
    sec = (lr-x0) * np.cos(0.5*arg1) * np.sin(arg2) / lr / np.cos(0.5*arg2) / arg2 
    #sec = (-1)**self.k * (lr-x0) * sin(self.k*pi-w)/(self.k*pi-w)/lr
    return (fir + f1*sec + 2*f2*np.cos(w)/lr) 
  
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


def reproduce_WenjunSun_PNAS():
  # read data to be fitted
  infile = open("WackWebb2.dat", 'r')
  h, k, q, F = read_data(infile, skip=1)
  h = np.array(h, int)
  k = np.array(k, int)
  q = np.array(q, float)
  F = np.array(F, float)

  # The following parameters reproduce one of Wenjun Sun's results
  # See RIPPLE~1/PROG_DIR/REFINE~1/REFINE.CMP
  mdf = MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=118, A=21.8, f1=1, f2=-9, rho_M=1, R_HM=2.1, 
            X_h=20.4, psi=0.1571)
  mdf.edp_par['f1'].vary = False
#  mdf.fit_lattice()
  mdf.fit_edp()
  mdf.report_edp()
  
  # The following parameters approximately reproduce one of Sun's results
  # for f1 and f2 free
  mdf = MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=103, A=20.0, f1=0.7, f2=-2, rho_M=1, R_HM=2.2, 
            X_h=20.4, psi=0.1571) 
  mdf.fit_edp()
  mdf.report_edp() 


if __name__ == "__main__":
#  reproduce_WenjunSun_PNAS()
  # read data to be fitted
  infile = open('ripple_066-067.dat', 'r')
  h, k, q, F = read_data(infile, skip=1)
  h = np.array(h, int)
  k = np.array(k, int)
  q = np.array(q, float)
  F = np.array(F, float) 
  
  # Work on SDF
  sdf = SDF(h, k, F, q, qx=None, qz=None, D=59.1, lambda_r=141.7, gamma=1.7174, 
            x0=123, A=22, rho_M=2, R_HM=4.2, X_h=16.1, psi=0.08727)   
  sdf.fit_lattice()
  sdf.fit_edp()
  sdf.report_edp()
  
  # Work on MDF
  mdf = MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=103, A=20.0, f1=1, f2=0, rho_M=1, R_HM=2.2, 
            X_h=20.4, psi=0.1571) 
#  mdf.fit_lattice()
#  mdf.fit_edp()
#  mdf.report_edp()  

  
  
  #Fourier_decomp(qx, qz, F, phase=None, N=201, xmin=-100, xmax=100, zmin=-100, zmax=100):
  
  # Optimization using the Ripple class
  #p = np.array([103, 18.6, 3, 20, 20.1, 0.0873])
  #obj.fit_SDF(par=p)
