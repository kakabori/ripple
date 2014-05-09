"""
reproduce_WenjunSun_PNAS
"""
from ripformfactor import *

# read data to be fitted
filename = 'form_factor/WackWebb2.dat'
h, k, q, F = read_data_4_columns(filename)
h = np.array(h, int)
k = np.array(k, int)
q = np.array(q, float)
F = np.array(F, float)

# The following parameters reproduce one of Wenjun Sun's results
# See RIPPLE~1/PROG_DIR/REFINE~1/REFINE.CMP
mdf = MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
          x0=118, A=21.8, f1=1, f2=-9, rho_M=1, R_HM=2.1, 
          X_h=20.4, psi=0.1571)
mdf.edp_par['f1'].vary = False
#mdf.fit_lattice()
mdf.fit_edp()
mdf.report_edp()
  
# The following parameters approximately reproduce one of Sun's results
# for f1 and f2 free
mdf = MDF(h, k, F, q, qx=None, qz=None, D=57.94, lambda_r=141.7, gamma=1.7174, 
          x0=103, A=20.0, f1=0.7, f2=-2, rho_M=1, R_HM=2.2, 
          X_h=20.4, psi=0.1571) 
mdf.fit_edp()
mdf.report_edp() 
