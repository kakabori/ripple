from ripintensity import *

###############################################################################
############################## __main__ #######################################
###############################################################################
if __name__ == "__main__":
  # read data to be fitted
  infilename = 'intensity/085_h5.dat'
  h, k, q, I, sigma, combined = read_data_5_columns(infilename)

###############################################################################
  # Work on SDF
  sdf = SDF(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714, 
            x0=103, A=18.6, 
            common_scale=51, R_HM=2.1, X_h=20.1, psi=0.08) 
  sdf.set_combined_peaks(combined)
#  sdf.set_mask(h=1, k=0, value=False)
#  sdf.set_mask(h=2, k=0, value=False)
#  sdf.set_mask(h=3, k=5, value=False)
#  sdf.set_mask(h=3, k=6, value=False)
#  sdf.set_mask(h=4, k=0, value=False)
#  sdf.fit_lattice()
  sdf.fit_edp()
  sdf.export_model_F("h5_F_sdf.dat")
#  sdf.report_edp()

###############################################################################
  # Work on MDF
  mdf = MDF(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714, 
            x0=104, A=22, f1=1, f2=-5, 
            common_scale=50, R_HM=2.2, X_h=20, psi=0.05) 
#  mdf.set_mask(h=1, k=0, value=False)
#  mdf.set_mask(h=2, k=0, value=False)
  mdf.fit_edp()
  mdf.export_model_F("h5_F_mdf.dat")
#  mdf.report_edp()  

###############################################################################
  # Work on S1G
  s1g = S1G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714, 
            x0=103, A=18.5, 
            rho_H1=10.77, Z_H1=20.86, sigma_H1=3.43,
            rho_M=9.23, sigma_M=1.67, psi=0.0873, common_scale=5)
#  s1g.set_combined_peaks(combined)
#  s1g.fit_lattice()
#  s1g.edp_par['rho_H1'].vary = False
#  s1g.edp_par['sigma_H1'].vary = False
  s1g.edp_par['rho_M'].vary = False
  s1g.edp_par['sigma_M'].vary = False 
  s1g.fit_edp()
  s1g.export_model_F("h5_F_s1g.dat")
#  s1g.report_edp()            

###############################################################################
  # Work on M1G
  m1g = M1G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
            x0=120, A=18.0, f1=0.6, f2=-1, 
            rho_H1=10.77, Z_H1=19.3, sigma_H1=3.43,
            rho_M=9.23, sigma_M=1.67, psi=0.3, common_scale=5)
#  m1g.fit_lattice()
  m1g.edp_par['rho_H1'].vary = False
  m1g.edp_par['sigma_H1'].vary = False
  m1g.edp_par['rho_M'].vary = False
  m1g.edp_par['sigma_M'].vary = False 
  m1g.fit_edp()
  m1g.export_model_F("h5_F_m1g.dat")
#  m1g.report_edp()   

###############################################################################
  # Work on S2G
  s2g = S2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
            x0=100, A=20, 
            rho_H1=9.91, Z_H1=19.45, sigma_H1=2.94,
            rho_H2=7.27, Z_H2=23.47, sigma_H2=1.47, 
            rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=5)
#  s2g.set_combined_peaks(combined)
#  s2g.fit_lattice()
  s2g.edp_par['rho_H1'].vary = False
  s2g.edp_par['sigma_H1'].vary = False
  s2g.edp_par['rho_H2'].vary = False
  s2g.edp_par['sigma_H2'].vary = False 
  s2g.edp_par['rho_M'].vary = False
  s2g.edp_par['sigma_M'].vary = False 
  s2g.fit_edp()
  s2g.export_model_F("h5_F_s2g.dat")
#  s2g.report_edp()

###############################################################################
  # Work on M2G
  m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
            x0=110, A=20, f1=0.6, f2=-1, 
            rho_H1=9.91, Z_H1=19.45, sigma_H1=2.94,
            rho_H2=7.27, Z_H2=23.47, sigma_H2=1.47, 
            rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=5)
#  m2g.set_combined_peaks(combined)
#  m2g.fit_lattice()
#  m2g.edp_par['rho_H1'].vary = False
  m2g.edp_par['sigma_H1'].vary = False
  m2g.edp_par['rho_H2'].vary = False
  m2g.edp_par['sigma_H2'].vary = False 
  m2g.edp_par['rho_M'].vary = False
  m2g.edp_par['sigma_M'].vary = False 
  m2g.fit_edp()
  m2g.export_model_F("h5_F_m2g.dat")
#  s2g.report_edp()
