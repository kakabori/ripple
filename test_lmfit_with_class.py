import numpy as np
import matplotlib as ppl
import scipy.optimize
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import lmfit
import math
  
class Exponential:
  def __init__(self, x, data):
    self.x = np.array(x, float)
    self.data = np.array(data, float)
    
    self.params = Parameters()
    self.params.add('A', value=1, vary=True)
    self.params.add('x0', value=1, vary=True)
    self.params.add('sigma', value=1, vary=True)

  
  def minimize(self):
    self.result = minimize(self.residuals, self.params)
    lmfit.report_fit(self.params)
    print("chisqr = {0:.3f}".format(self.result.chisqr))

    
  def residuals(self, params):
    A = params['A'].value
    x0 = params['x0'].value
    s = params['sigma'].value 
  
    model = self.model(A, x0, s)
    return (model - self.data)

  
  def model(self, A=None, x0=None, s=None, x=None):
    if x is None:
      x = self.x
    if s is None:
      s = self.params['sigma'].value
    if x0 is None:
      x0 = self.params['x0'].value
    if A is None:
      A = self.params['A'].value  
    
    return A * np.exp(-(x-x0)**2/2/s)
    
    
x = np.linspace(-8, 5, 51)
y = np.array(100*np.exp(-(x+4)**2/2/3.352))
tmp = Exponential(x, y)
tmp.minimize()


