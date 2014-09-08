from ripplefit import *

# read data to be fitted
infilename = 'intensity/085_h9_ver4.dat'
h, k, q, I, sigma, combined = read_data_5_columns(infilename)

###############################################################################
# Work on M2G
m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
          xM=90, A=25, f1=1.5, f2=-20, 
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

m2g.report_edp()
#m2g.export_model_F("fits/fit5_F.txt")
#m2g.export_model_I("fits/fit5_I.txt")
#m2g.export_2D_edp("fits/fit5_2D_edp.txt")
#m2g.export_params("fits/fit5_params.txt")
#m2g.export_EDP("fits/fit5_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_1.txt", center=(72.5,0), angle=-15, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_2.txt", center=(72.5,0), angle=-10, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_3.txt", center=(72.5,0), angle=-8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_4.txt", center=(72.5,0), angle=-6, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_5.txt", center=(72.5,0), angle=-4, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_6.txt", center=(72.5,0), angle=-2, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_1D_minor_7.txt", center=(72.5,0), angle=0, length=100, stepsize=0.1)
m2g.export_headgroup_positions("fits/fit5_headgroup.txt")
m2g.export_methyl_positions("fits/fit5_methyl.txt")
m2g.export_phases("fits/fit5_phases.txt")
