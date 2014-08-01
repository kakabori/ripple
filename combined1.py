""" 
EDP created based on the phases from Fit5 but with some orders reversed
"""

from ripintensity import *
import copy

# read data to be fitted
infilename = 'intensity/combined.dat'
h, k, q, I, sigma, combined = read_data_5_columns(infilename)

tmp = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
          x0=98.91, A=20.6, f1=0.575, f2=-2.92, 
          rho_H1=5.385, Z_H1=19.75, sigma_H1=3.43,
          rho_H2=5.385, Z_H2=19.75, sigma_H2=3.43, 
          rho_M=9.23, sigma_M=1.67, psi=0.261, common_scale=1.02)
tmp.edp_par['f1'].vary = True
tmp.edp_par['f2'].vary = True
tmp.edp_par['Z_H2'].vary = True
tmp.edp_par['rho_H1'].vary = False
tmp.edp_par['sigma_H1'].vary = False
tmp.edp_par['rho_H2'].vary = False
tmp.edp_par['sigma_H2'].vary = False
tmp.edp_par['rho_M'].vary = False
tmp.edp_par['sigma_M'].vary = False
tmp.fit_edp()
  
tmp.edp_par['sigma_H1'].vary = True
tmp.edp_par['sigma_H2'].vary = True
tmp.edp_par['sigma_M'].vary = True
tmp.fit_edp()

tmp.edp_par['rho_H1'].vary = True
tmp.edp_par['rho_H2'].vary = True
tmp.edp_par['rho_M'].vary = True
tmp.fit_edp()

###############################################################################  
# read data to be fitted
infilename = 'intensity/085_h9_ver4.dat'
h, k, q, I, sigma, combined = read_data_5_columns(infilename)

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

# Copy phases of m2g to tmp
tmp.phase = m2g.phase

# Now flip some of the phases
tmp.phase[(tmp.h==3)&(tmp.k==0)] *= -1
tmp.phase[(tmp.h==6)&(tmp.k==0)] *= -1
tmp.phase[(tmp.h==6)&(tmp.k==1)] *= -1
tmp.phase[(tmp.h==6)&(tmp.k==3)] *= -1
tmp.phase[(tmp.h==6)&(tmp.k==4)] *= -1
tmp.phase[(tmp.h==9)&(tmp.k==0)] *= -1

tmp.export_phases("fits/combined1_phases.txt")
tmp.export_model_F("fits/combined1_F.txt")
tmp.export_model_I("fits/combined1_I.txt")
tmp.export_2D_edp("fits/combined1_2D_edp.txt")
tmp.export_params("fits/combined1_params.txt")
tmp.export_angle("fits/combined1_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
tmp.export_angle("fits/combined1_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
tmp.export_headgroup_positions("fits/combined1_headgroup.txt")
