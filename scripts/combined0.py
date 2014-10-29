"""
This script fits data where lower orders are replaced with Wack and Webb
data set.
"""

from ripintensity import *

# read data to be fitted
infilename = 'intensity/combined.dat'
h, k, q, I, sigma, combined = read_data_5_columns(infilename)

###############################################################################
# Work on M2G
m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
          x0=98.91, A=20.6, f1=0.575, f2=-2.92, 
          rho_H1=5.385, Z_H1=19.75, sigma_H1=3.43,
          rho_H2=5.385, Z_H2=19.75, sigma_H2=3.43, 
          rho_M=9.23, sigma_M=1.67, psi=0.261, common_scale=1.02)
m2g.edp_par['f1'].vary = True
m2g.edp_par['f2'].vary = True
m2g.edp_par['Z_H2'].vary = True
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
  
m2g.report_edp()
m2g.export_phases("fits/combined0_phases.txt")
m2g.export_model_F("fits/combined0_F.txt")
m2g.export_model_I("fits/combined0_I.txt")
m2g.export_2D_edp("fits/combined0_2D_edp.txt")
m2g.export_params("fits/combined0_params.txt")
m2g.export_angle("fits/combined0_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_angle("fits/combined0_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.export_headgroup_positions("fits/combined0_headgroup.txt")