"""
This script creates and outputs a bunch of EDP's with various combinations
of the suspected phases with base phases given by those from Fit5.
"""

from ripintensity import *

# read data to be fitted
infilename = 'intensity/085_h9_ver4.dat'
h, k, q, I, sigma, combined = read_data_5_columns(infilename)

###############################################################################
# Work on M2G
m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
          x0=90, A=25, f1=1.5, f2=-20, 
          rho_H1=9.91, Z_H1=20, sigma_H1=2.94,
          rho_H2=7.27, Z_H2=20, sigma_H2=1.47, 
          rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=3)
#m2g.set_combined_peaks(combined)
#m2g.fit_lattice()
m2g.edp_par['f2'].vary = True
m2g.edp_par['rho_H1'].vary = False
m2g.edp_par['sigma_H1'].vary = False
m2g.edp_par['rho_H2'].vary = False
m2g.edp_par['sigma_H2'].vary = False
m2g.edp_par['rho_M'].vary = False
m2g.edp_par['sigma_M'].vary = False
m2g.fit_edp()
  
m2g.edp_par['sigma_H1'].vary = True
m2g.edp_par['sigma_H2'].vary = True
m2g.edp_par['sigma_M'].vary = True
m2g.fit_edp()
  
m2g.edp_par['rho_H1'].vary = True
m2g.edp_par['rho_H2'].vary = True
m2g.edp_par['rho_M'].vary = True
m2g.fit_edp()

m2g.phase[(m2g.h==1)&(m2g.k==3)] *= -1
m2g.phase[(m2g.h==3)&(m2g.k==0)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==0)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==1)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==3)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==4)] *= -1
m2g.phase[(m2g.h==9)&(m2g.k==0)] *= -1
#m2g.export_2D_edp("fits/fit22_2D_edp.txt")
#m2g.export_angle("fits/fit22_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_angle("fits/fit22_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.export_headgroup_positions("fits/fit22_headgroup.txt")
m2g.export_phases("fits/fit22_phases.txt")
m2g.phase[(m2g.h==1)&(m2g.k==3)] *= -1
m2g.phase[(m2g.h==3)&(m2g.k==0)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==0)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==1)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==3)] *= -1
m2g.phase[(m2g.h==6)&(m2g.k==4)] *= -1
m2g.phase[(m2g.h==9)&(m2g.k==0)] *= -1
