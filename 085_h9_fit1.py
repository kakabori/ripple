from ripintensity import *

###############################################################################
############################## __main__ #######################################
###############################################################################
if __name__ == "__main__":
  # read data to be fitted
  infilename = 'intensity/085_h9_ver4.dat'
  h, k, q, I, sigma, combined = read_data_5_columns(infilename)

###############################################################################
  # Work on M1G
  m1g = M1G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
            x0=110.5, A=25.78, f1=0.8, f2=0, 
            rho_H_major=10.77, rho_H_minor=10.77,
            Z_H_major=20.08, Z_H_minor=20.08,
            sigma_H_major=3.43, sigma_H_minor=3.43,
            rho_M_major=9.23, rho_M_minor=9.23,
            sigma_M_major=1.67, sigma_M_minor=1.67,
            psi_major=0.32, psi_minor=0.32,
            common_scale=3)
#  m1g.fit_lattice()
  m1g.edp_par['f1'].vary = True
  m1g.edp_par['f2'].vary = False
  m1g.edp_par['rho_H_major'].vary = False
  m1g.edp_par['rho_H_minor'].vary = False
  m1g.link_rho_H = True
  m1g.edp_par['Z_H_major'].vary =True
  m1g.edp_par['Z_H_minor'].vary = False
  m1g.link_Z_H = True
  m1g.edp_par['sigma_H_major'].vary = False
  m1g.edp_par['sigma_H_minor'].vary = False
  m1g.link_sigma_H = True
  m1g.edp_par['rho_M_major'].vary = False
  m1g.edp_par['rho_M_minor'].vary = False
  m1g.link_rho_M = True
  m1g.edp_par['sigma_M_major'].vary = False
  m1g.edp_par['sigma_M_minor'].vary = False
  m1g.link_sigma_M = True
  m1g.edp_par['psi_major'].vary = True
  m1g.edp_par['psi_minor'].vary = False
  m1g.link_psi = True
  m1g.fit_edp()
  m1g.export_model_F("fit1_F.txt")
  m1g.export_model_I("fit1_I.txt")
  m1g.export_2D_edp("fit1_2D_edp.txt")
  m1g.export_params("fit1_params.txt")
  m1g.export_angle("fit1_1D_major.txt", center=(0,0), angle=-11.8, length=80, stepsize=0.1)
  m1g.export_angle("fit1_1D_minor.txt", center=(72.5,0), angle=27.1, length=80, stepsize=0.1)
  m1g.report_edp()   

###############################################################################
  # Work on M2G
  m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
            x0=90, A=25, f1=0.5, f2=-10, 
            rho_H1=9.91, Z_H1=15.45, sigma_H1=2.94,
            rho_H2=7.27, Z_H2=23.47, sigma_H2=1.47, 
            rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=3)
  m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
            x0=120, A=18, f1=0.6, f2=-1, 
            rho_H1=5.385, Z_H1=19.3, sigma_H1=3.43,
            rho_H2=5.385, Z_H2=19.3, sigma_H2=3.43, 
            rho_M=9.23, sigma_M=1.67, psi=0.3, common_scale=1)
#  m2g.set_combined_peaks(combined)
#  m2g.fit_lattice()
  m2g.edp_par['rho_H1'].vary = False
  m2g.edp_par['sigma_H1'].vary = False
  m2g.edp_par['rho_H2'].vary = False
  m2g.edp_par['sigma_H2'].vary = False
  m2g.edp_par['rho_M'].vary = False
  m2g.edp_par['sigma_M'].vary = False
#  m2g.fit_edp()
#  m2g.export_model_F("h9_F_m2g_sH2_free.dat")
#  m2g.export_model_I("h9_I_m2g_sH2_free.dat")
#  m2g.export_2D_edp("h9_2D_edp_m2g_sH2_free.dat")
#  m2g.report_edp()

###############################################################################
  # Work on S2G
#  s2g = S2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
#            x0=100, A=25, 
#            rho_H1=9.91, Z_H1=19.45, sigma_H1=2.94,
#            rho_H2=7.27, Z_H2=23.47, sigma_H2=1.47, 
#            rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=5)
#  s2g.set_combined_peaks(combined)
#  s2g.fit_lattice()
#  s2g.edp_par['rho_H1'].vary = False
#  s2g.edp_par['sigma_H1'].vary = False
#  s2g.edp_par['rho_H2'].vary = False
#  s2g.edp_par['sigma_H2'].vary = False 
#  s2g.edp_par['rho_M'].vary = False
#  s2g.edp_par['sigma_M'].vary = False 
#  s2g.fit_edp()
#  s2g.export_model_F("h9_F_s2g.dat")
#  s2g.report_edp()

###############################################################################
  # Work on S1G
#  s1g = S1G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714, 
#            x0=103, A=18.5, 
#            rho_H1=10.77, Z_H1=20.86, sigma_H1=3.43,
#            rho_M_major=9.23, sigma_M=1.67, psi=0.0873, common_scale=5,
#            rho_M_minor=9.23)
#  s1g.set_combined_peaks(combined)
#  s1g.fit_lattice()
#  s1g.edp_par['rho_H1'].vary = False
#  s1g.edp_par['sigma_H1'].vary = False
#  s1g.edp_par['rho_M_major'].vary = False
#  s1g.edp_par['rho_M_minor'].vary = False
#  s1g.edp_par['sigma_M'].vary = False 
#  s1g.fit_edp()
#  s1g.export_model_F("h9_F_s1g.dat")
#  s1g.report_edp()  

###############################################################################
  # Work on SDF
#  sdf = SDF(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714, 
#            x0=103, A=18.6, 
#            common_scale=51, R_HM=2.1, X_h=20.1, psi=0.08) 
#  sdf.set_combined_peaks(combined)
#  sdf.set_mask(h=1, k=0, value=False)
#  sdf.set_mask(h=2, k=0, value=False)
#  sdf.set_mask(h=3, k=5, value=False)
#  sdf.set_mask(h=3, k=6, value=False)
#  sdf.set_mask(h=4, k=0, value=False)
#  sdf.fit_lattice()
#  sdf.fit_edp()
#  sdf.export_model_F("h7_F_sdf.dat")
#  sdf.report_edp()

###############################################################################
  # Work on MDF
#  mdf = MDF(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714, 
#            x0=65.5, A=22.9, f1=1, f2=0, 
#            common_scale=70, R_HM=0.8, X_h=20.7, psi=0.4) 
#  mdf.set_mask(h=1, k=0, value=False)
#  mdf.set_mask(h=2, k=0, value=False)
#  mdf.edp_par['f1'].vary = False
#  mdf.edp_par['f2'].vary = True
#  mdf.edp_par['common_scale'].max = 100
#  mdf.fit_edp()
#  mdf.export_model_F("h7_F_mdf.dat")
#  mdf.report_edp()       
