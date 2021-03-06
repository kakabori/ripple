import sys
import numpy as np
from numpy import pi, sin, cos, tan, exp, sqrt
#import matplotlib as ppl
import scipy.optimize
#import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit
from scipy import ndimage
import math
import random

# A module-level global variable
wavelength = 1.175

def read_data_3_columns(filename):
    # Process comment and header lines
    fileobj = open(filename, 'r')
    comment_and_header(fileobj)  
    # Go through data points  
    hl = []; kl = []; pl = []; 
    lines = fileobj.readlines()
    for line in lines:
        # This ignores an empty line
        line = line.rstrip()
        if not line: 
            continue
        h, k, p = line.split()
        h = int(h)
        k = int(k)
        p = int(p)
        hl.append(h); kl.append(k); pl.append(p);       
    return hl, kl, pl
    
def read_data_4_columns(filename):
    """Read a four-column ASCII file and parse each column into a python list.
    Lines starting with # will be ignored, i.e., # signals a comment line.
    
    filename: input file name
    
    The input file must be formatted as "h k F sigma_F", where 
    h, k : ripple main and side peak index
    F : form factor
    sigma_F : uncertainty in F
    """
    # Process comment and header lines
    fileobj = open(filename, 'r')
    comment_and_header(fileobj)      
    # Go through data points  
    hl = []; kl = []; Fl = []; sl =[]; 
    lines = fileobj.readlines()
    for line in lines:
        # This ignores an empty line
        line = line.rstrip()
        if not line: 
            continue
        h, k, F, s = line.split()
        h = int(h)
        k = int(k)
        F = float(F)
        s = float(s)
        hl.append(h); kl.append(k); Fl.append(F); sl.append(s)      
    return hl, kl, Fl, sl
     
def read_data_5_columns(filename):
    """Read a five-column ASCII file and parse each column into a python list.
    Lines starting with # will be ignored, i.e., # signals a comment line.
    
    filename: input file name
    
    The input file must be formatted as "h k q I sigma", where 
    h, k: ripple main and side peak index
    q: magnitude of scattering vector, q
    I: observed intensity
    sigma: uncertainty in intensity
    
    Example
    =======
    # Comment goes here
    # Another comment line
    h  k      q      I  sigma
    1 -1  0.107   78.9    8.1
    1  0  0.100  100.0   10.0
    2  0  0.200   45.6    6.9
    2  1  0.205   56.7    8.0
    """
    # Process comment and header lines
    fileobj = open(filename, 'r')
    comment_and_header(fileobj)   
    # Go through data points  
    hl = []; kl = []; ql = []; Il = []; sl =[]
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
        q = float(q)
        hl.extend(h); kl.extend(k); ql.append(q); Il.append(I); sl.append(s)     
    return hl, kl, ql, Il, sl

def comment_and_header(fileobj):
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
    """
    def __init__(self, h, k, I, sigma_I, D, lambda_r, gamma, F=None, sigma_F=None):
        self.h = np.array(h, int)
        self.k = np.array(k, int)
        self.I = np.array(I, float)
        self.sigma = np.array(sigma_I, float)
        self.latt_par = Parameters()
        self.latt_par.add('D', value=D, vary=True)
        self.latt_par.add('lambda_r', value=lambda_r, vary=True)
        self.latt_par.add('gamma', value=gamma, vary=True)
        self._set_qxqz()
        self.F = np.array(F, float)
        self.sigma_F = np.array(sigma_F, float)
    
    def update(self, h, k, I, sigma_I):
        self.h = np.array(h, int)
        self.k = np.array(k, int)
        self.I = np.array(I, float)
        self.sigma = np.array(sigma_I, float)
               
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
        exp_F = sqrt(self.I)
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
        """Start a non-linear least squared fit for electron density profile.
        After the fit is completed, the form factor is updated."""
        self._set_qxqz()
        self.edp = minimize(self._residual_edp, self.edp_par)
        phase = np.sign(self._model_F())
        phase = phase.astype(int)
        self.F = phase * sqrt(self.I)
        expF_10 = np.abs(self.F[(self.h==1)&(self.k==0)])
        self.F = self.F / expF_10 * 100
    
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
        exp_F = sqrt(self.I)
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
        """Export an ASCII file containing phase factors predicted
        by a model"""
        phase = np.sign(self._model_F())
        phase = phase.astype(int)
        with open(filename, 'w') as f:
            f.write("h k phase\n")
            for a, b, c in zip(self.h, self.k, phase):
                f.write("{0: 1d} {1: 1d} {2: 1d}\n".format(a, b, c))

    def get_phase(self, h, k):
        """Get the phase factor of the (h,k) order"""
        return self.phase[(self.h==h)&(self.k==k)]

    def get_normalized_exp_F(self):
        exp_F = sqrt(self.I)
        # exp_F is normalized at (h=1,k=0) peak
        expF_10 = exp_F[(self.h==1)&(self.k==0)]
        exp_F = exp_F / expF_10 * 100
        return exp_F
           
                                
###############################################################################
class Sawtooth(BaseRipple):
  def __init__(self, h, k, I, sigma, D=57.8, lambda_r=145, gamma=1.71, 
               xM=100, A=20):
    super(Sawtooth, self).__init__(h, k, I, sigma, D, lambda_r, gamma)
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
  def __init__(self, h, k, I, sigma, D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20.27, 
               rho_H1=2.21, Z_H1=20.00, sigma_H1=3.33,
               rho_H2=2.22, Z_H2=22.22, sigma_H2=3.33,
               rho_M=1, sigma_M=3, psi=0.087, common_scale=0.1):
    super(S2G, self).__init__(h, k, I, sigma, D, lambda_r, gamma, xM, A)
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
  def __init__(self, h=[], k=[], I=[], sigma=[],
               D=58, lambda_r=140, gamma=1.7, 
               xM=100, A=20, f1=1, f2=0, 
               rho_H1=2.21, Z_H1=20.24, sigma_H1=3.33,
               rho_H2=2.22, Z_H2=20.22, sigma_H2=3.33, 
               rho_M=1, sigma_M=3, psi=0.087, common_scale=0.1):
    super(M2G, self).__init__(h, k, I, sigma, D, lambda_r, gamma, xM, A, 
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
   

