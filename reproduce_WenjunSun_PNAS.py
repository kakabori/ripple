"""
reproduce_WenjunSun_PNAS
"""
import ripformfactor as rf

# read data to be fitted
filename = 'form_factor/WackWebb2.dat'
h, k, q, F = rf.read_data_4_columns(filename)
h = np.array(h, int)
k = np.array(k, int)
q = np.array(q, float)
F = np.array(F, float)

# The following parameters reproduce one of Wenjun Sun's results
# See RIPPLE~1/PROG_DIR/REFINE~1/REFINE.CMP
mdf = rf.MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=118, A=21.8, f1=1, f2=-9, rho_M=1, R_HM=2.1, 
            X_h=20.4, psi=0.1571)
mdf.edp_par['f1'].vary = False
#mdf.fit_lattice()
mdf.fit_edp()
mdf.report_edp()
  
# The following parameters approximately reproduce one of Sun's results
# for f1 and f2 free
mdf = rf.MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=103, A=20.0, f1=0.7, f2=-2, rho_M=1, R_HM=2.2, 
            X_h=20.4, psi=0.1571) 
mdf.fit_edp()
mdf.report_edp() 

import ripintensity as ri
# read data again, but in intensity
filename = 'intensity/WackWebb2.dat'
h, k, q, I, sigma = read_data_5_columns(filename)
h = np.array(h, int)
k = np.array(k, int)
q = np.array(q, float)
I = np.array(I, float)
sigma = np.array(sigma, float)

# Work on S1G
s1g = S1G(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.7174,
          x0=100, A=20, rho_H1=10.77, Z_H1=20.86, sigma_H1=3.43,
          rho_M=9.23, sigma_M=1.67, psi=0.1, DeltaRho=1)
s2g.fit_lattice()
s2g.edp_par['rho_H1'].vary = False
s2g.edp_par['sigma_H1'].vary = False
s2g.edp_par['rho_M'].vary = False
s2g.edp_par['sigma_M'].vary = False 
s2g.fit_edp()
s2g.report_edp()            
            
# Work on S2G
s2g = S2G(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.7174,
          x0=100, A=20, rho_H1=9.91, Z_H1=19.45, sigma_H1=2.94,
          rho_H2=7.27, Z_H2=23.47, sigma_H2=1.47, 
          rho_M=10.91, sigma_M=1.83, psi=0.1, DeltaRho=1)
s2g.fit_lattice()
s2g.edp_par['rho_H1'].vary = False
s2g.edp_par['sigma_H1'].vary = False
s2g.edp_par['rho_H2'].vary = False
s2g.edp_par['sigma_H2'].vary = False 
s2g.edp_par['rho_M'].vary = False
s2g.edp_par['sigma_M'].vary = False 
s2g.fit_edp()
s2g.report_edp()
