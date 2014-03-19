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
    """
    params is a dummy variable
    """
    model = self.model()
    data = self.data
    return (model - data)
  
  def model(self):
    A = self.params['A'].value
    x0 = self.params['x0'].value
    s = self.params['sigma'].value    
    return A * np.exp(-(self.x-x0)**2/2/s)
    
    
x = np.linspace(-8, 5, 51)
y = np.array(120*np.exp(-(x+4)**2/2/3.352))
tmp = Exponential(x, y)
tmp.minimize()


