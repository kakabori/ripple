from ripplefit import *

# read data to be fitted
infilename = 'intensity/085_h9_ver6.dat'
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

m2g.edp_par['rho_H1'].vary = True
m2g.edp_par['rho_H2'].vary = True
m2g.edp_par['rho_M'].vary = True
m2g.fit_edp()

m2g.edp_par['sigma_H1'].vary = True
m2g.edp_par['sigma_H2'].vary = True
m2g.edp_par['sigma_M'].vary = True
m2g.fit_edp()
  
#m2g.report_edp()
#m2g.export_model_F("fits/fit7_F.txt")
#m2g.export_model_I("fits/fit7_I.txt")
#m2g.export_2D_edp("fits/fit7_2D_edp.txt")
#m2g.export_params("fits/fit7_params.txt")
#m2g.export_angle("fits/fit7_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_angle("fits/fit7_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
#m2g.export_headgroup_positions("fits/fit7_headgroup.txt")
#m2g.export_methyl_positions("fits/fit7_methyl.txt")
#m2g.export_phases("fits/fit7_phases.txt")

m2g.export_EDP("fits/fit7_EDP_0.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit7_EDP_1.txt", center=(10,2.1), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit7_EDP_2.txt", center=(20,4.3), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit7_EDP_3.txt", center=(30,6.4), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit7_EDP_4.txt", center=(40,8.5), angle=-11.8, length=100, stepsize=0.1)

m2g.export_EDP_between_two_points("fits/fit7_EDP_11.txt", start=(-40,-1), end=(40,16), N=161)
m2g.export_EDP_between_two_points("fits/fit7_EDP_12.txt", start=(-40,21), end=(40,38), N=161)

# Plot through water region
m2g.export_EDP_between_two_points("fits/PF7_EDP_13.txt", start=(-145,29), end=(-95,39), N=51)
m2g.export_EDP_between_two_points("fits/PF7_EDP_14.txt", start=(-95,39), end=(-55,18), N=41)
m2g.export_EDP_between_two_points("fits/PF7_EDP_15.txt", start=(-55,18), end=(50,39), N=106)
m2g.export_EDP_between_two_points("fits/PF7_EDP_16.txt", start=(50,39), end=(90,18), N=41)
m2g.export_EDP_between_two_points("fits/PF7_EDP_17.txt", start=(90,18), end=(145,29), N=56)
