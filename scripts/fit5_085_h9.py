import sys
sys.path.append('/home/unkokusei/prog/ripple/src')

from ripple import *

# read data to be fitted
#infilename = 'intensity/085_h9_ver6.dat'
infilename = '../intensity/085_h9_ver4.dat'
h, k, q, I, sigma = read_data_5_columns(infilename)

###############################################################################
# Work on M2G
m2g = M2G(h, k, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
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

#m2g.report_edp()
#m2g.export_model_F("fits/fit5_F.txt")
#m2g.export_model_I("fits/fit5_I.txt")
#m2g.export_2D_edp("fits/fit5_2D_edp.txt")
#m2g.export_params("fits/fit5_params.txt")
#m2g.export_angle("fits/fit5_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_angle("fits/fit5_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
#m2g.export_headgroup_positions("fits/fit5_headgroup.txt")
#m2g.export_methyl_positions("fits/fit5_methyl.txt")
#m2g.export_phases("fits/fit5_phases.txt")

# Plot EDP on a line with xiM = 12 degrees
#m2g.export_EDP("fits/fit5_EDP_0.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_1.txt", center=(10,2.1), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_2.txt", center=(20,4.3), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_3.txt", center=(30,6.4), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_4.txt", center=(40,8.5), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_5.txt", center=(-10,-2.1), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_6.txt", center=(-20,-4.3), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_7.txt", center=(-30,-6.4), angle=-11.8, length=100, stepsize=0.1)
#m2g.export_EDP("fits/fit5_EDP_8.txt", center=(-40,-8.5), angle=-11.8, length=100, stepsize=0.1)

# Plot through hydrocarbon region
#m2g.export_EDP_between_two_points("fits/fit5_EDP_11.txt", start=(-40,-1), end=(40,16), N=161)

# Plot through water region
#m2g.export_EDP_between_two_points("fits/fit5_EDP_13.txt", start=(-145,29), end=(-95,39), N=51)
#m2g.export_EDP_between_two_points("fits/fit5_EDP_14.txt", start=(-95,39), end=(-55,18), N=41)
#m2g.export_EDP_between_two_points("fits/fit5_EDP_15.txt", start=(-55,18), end=(50,39), N=106)
#m2g.export_EDP_between_two_points("fits/fit5_EDP_16.txt", start=(50,39), end=(90,18), N=41)
#m2g.export_EDP_between_two_points("fits/fit5_EDP_17.txt", start=(90,18), end=(145,29), N=56)

# Apply mask
m2g.export_EDP("fits/fit5_EDP_20.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_EDP_30.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.edm.data.set_mask(True, h=9)
m2g.export_EDP("fits/fit5_EDP_21.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_EDP_31.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.edm.data.set_mask(True, h=7)
m2g.export_EDP("fits/fit5_EDP_22.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_EDP_32.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.edm.data.set_mask(True, h=6)
m2g.export_EDP("fits/fit5_EDP_23.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_EDP_33.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.edm.data.set_mask(True, h=5)
m2g.export_EDP("fits/fit5_EDP_24.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_EDP_34.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
m2g.edm.data.set_mask(True, h=4)
m2g.export_EDP("fits/fit5_EDP_25.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
m2g.export_EDP("fits/fit5_EDP_35.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
