import numpy as np
import matplotlib.pyplot as plt

def strip(x):
  global rho_w, rho_CH2, x0, sigma_H
  w = 2.36*sigma_H
  ret = np.zeros(len(x))
  ret[(x>=0) & (x<=x0-w/2)] = rho_CH2
  ret
  [(x>=x0-w/2) & (x<=x0+w/2)] = f(x[(x>=x0-w/2) & (x<=x0+w/2)], w) 
  ret[x>=x0+w/2] = rho_w

  return ret    
    
    
def _f(x, w):
  global rho_w, rho_CH2, x0
  return (rho_w-rho_CH2)/2 * np.sin(np.pi*(x-x0)/w) + (rho_w+rho_CH2)/2
  

def gauss_part(x):
  global sigma_H, x0, A_H
  global sigma_M, A_M
  global rho_w
  head = A_HG * np.exp(-(x-x0)**2/2/sigma_HG/sigma_HG)
  methyl = A_M * np.exp(-x**2/2/sigma_M/sigma_M)
  return (head + methyl)
  

rho_w = 0.333
rho_CH2 = 0.222
x0 = 18
sigma_H = 3
A_H = 0.5
sigma_M = 2
A_M = -0.1
    
if __name__ == "__main__":
  x = np.linspace(0, 35, 100)
  y = gauss_part(x) + strip(x)
  plt.plot(x, y)
  plt.show()
