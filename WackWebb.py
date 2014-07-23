from ripintensity import *

###############################################################################
############################## __main__ #######################################
###############################################################################
if __name__ == "__main__":
  # read data to be fitted
  infilename = 'intensity/WackWebb.dat'
  h, k, q, I, sigma, combined = read_data_5_columns(infilename)

###############################################################################
  # Work on M1G
#  m1g = M1G(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.717,
#            x0=115, A=24, f1=0.7, f2=-10, 
#            rho_H_major=10.77, rho_H_minor=10.77,
#            Z_H_major=21, Z_H_minor=21,
#            sigma_H_major=3.43, sigma_H_minor=3.43,
#            rho_M_major=9.23, rho_M_minor=9.23,
#            sigma_M_major=1.67, sigma_M_minor=1.67,
#            psi_major=0.3, psi_minor=0.3,
#            common_scale=3)
  m1g = M1G(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.717,
            x0=103, A=19, f1=0.6, f2=-1, 
            rho_H_major=10.77, rho_H_minor=10.77,
            Z_H_major=20.08, Z_H_minor=20.08,
            sigma_H_major=3.43, sigma_H_minor=3.43,
            rho_M_major=9.23, rho_M_minor=9.23,
            sigma_M_major=1.67, sigma_M_minor=1.67,
            psi_major=0.157, psi_minor=0.157,
            common_scale=3)
#  m1g.fit_lattice()
  m1g.edp_par['f1'].vary = True
  m1g.edp_par['f2'].vary = True
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
  m1g.export_model_F("fit9_F.txt")
  m1g.export_model_I("fit9_I.txt")
  m1g.export_2D_edp("fit9_2D_edp.txt")
  m1g.export_params("fit9_params.txt")
  m1g.export_angle("fit9_1D_major.txt", center=(0,0), angle=-11.8, length=100, stepsize=0.1)
  m1g.export_angle("fit9_1D_minor.txt", center=(72.5,0), angle=27.1, length=100, stepsize=0.1)
  m1g.report_edp()   
