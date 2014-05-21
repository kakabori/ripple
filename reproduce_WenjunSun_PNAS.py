"""
reproduce_WenjunSun_PNAS
"""
import ripformfactor as rf
import ripintensity as ri
import numpy as np

# read data to be fitted
filename = 'form_factor/WackWebb2.dat'
h, k, q, F = rf.read_data_4_columns(filename)
h = np.array(h, int)
k = np.array(k, int)
q = np.array(q, float)
F = np.array(F, float)

###############################################################################
sdf = rf.SDF(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174, 
             x0=102.8, A=21.1, common_scale=1.04, R_HM=2.52, 
             X_h=19.7, psi=0.0823) 
sdf.set_qxqz()
#sdf.fit_edp()
#print("Report on SDF form factor")
#sdf.report_edp()

###############################################################################  
mdf = rf.MDF(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=103, A=20.0, f1=0.7, f2=-2, common_scale=1, R_HM=2.2, 
            X_h=20.4, psi=0.1571) 
mdf.set_qxqz()
#mdf.fit_edp()
#print("Report on MDF form factor")
#mdf.report_edp() 

###############################################################################
# The following parameters reproduce one of Wenjun Sun's results
# See RIPPLE~1/PROG_DIR/REFINE~1/REFINE.CMP
mdf = rf.MDF(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=118, A=21.8, f1=1, f2=-9, common_scale=1, R_HM=2.1, 
            X_h=20.4, psi=0.1571)
mdf.edp_par['f1'].vary = False
mdf.set_qxqz()
#mdf.fit_edp()
#mdf.report_edp()

###############################################################################  
# The following parameters approximately reproduce one of Sun's results
# for f1 and f2 free
mdf = rf.MDF(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174, 
            x0=103, A=20.0, f1=0.7, f2=-2, common_scale=1, R_HM=2.2, 
            X_h=20.4, psi=0.1571) 
mdf.edp_par['x0'].vary = False
mdf.edp_par['A'].vary = False
mdf.edp_par['f1'].vary = False
mdf.edp_par['f2'].vary = False
mdf.edp_par['common_scale'].vary = True
mdf.edp_par['R_HM'].vary = True
mdf.edp_par['X_h'].vary = False
mdf.edp_par['psi'].vary = False
mdf.set_qxqz()
#mdf.fit_edp()
#mdf.report_edp() 

###############################################################################
s1g = rf.S1G(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174,
             x0=103, A=18.5, 
             rho_H1=10.77, Z_H1=20.86, sigma_H1=3.43,
             rho_M=9.23, sigma_M=1.67, psi=0.0873, common_scale=1)
s1g.edp_par['rho_H1'].vary = False
s1g.edp_par['sigma_H1'].vary = False
s1g.edp_par['rho_M'].vary = False
s1g.edp_par['sigma_M'].vary = False 
s1g.set_qxqz()
#print("Report on S1G form factor")
#s1g.fit_edp()
#s1g.report_edp()            
            
###############################################################################
# Work on M1G
m1g = rf.M1G(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174,
             x0=103, A=19.0, f1=0.6, f2=-1, 
             rho_H1=10.77, Z_H1=19.3, sigma_H1=3.43,
             rho_M=9.23, sigma_M=1.67, psi=0.157, common_scale=1)
m1g.edp_par['rho_H1'].vary = False
m1g.edp_par['sigma_H1'].vary = False
m1g.edp_par['rho_M'].vary = False
m1g.edp_par['sigma_M'].vary = False 
m1g.set_qxqz()
print("Report on M1G form factor")
m1g.fit_edp()
m1g.report_edp()   
            
###############################################################################            
# Work on S2G. This model was not used in WenjunSun's PNAS
s2g = rf.S2G(h, k, q, F, D=57.94, lambda_r=141.7, gamma=1.7174,
             x0=100, A=20, rho_H1=9.91, Z_H1=19.45, sigma_H1=2.94,
             rho_H2=7.27, Z_H2=23.47, sigma_H2=1.47, 
             rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=1)
s2g.edp_par['rho_H1'].vary = False
s2g.edp_par['sigma_H1'].vary = False
s2g.edp_par['rho_H2'].vary = False
s2g.edp_par['sigma_H2'].vary = False 
s2g.edp_par['rho_M'].vary = False
s2g.edp_par['sigma_M'].vary = False 
#s2g.fit_edp()
#s2g.report_edp()

###############################################################################
############ Check that ripintensity.py works #################################
###############################################################################
# read data to be fitted
infilename = 'intensity/WackWebb.dat'
h, k, q, I, sigma, combined = ri.read_data_5_columns(infilename)

###############################################################################
# Work on SDF
sdf2 = ri.SDF(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.7174, 
              x0=102.8, A=21.1, common_scale=1.04, 
              R_HM=2.52, X_h=19.7, psi=0.0823) 
#sdf2.fit_edp()
#print("Report on SDF intensity")
#sdf2.report_edp()

###############################################################################
# Work on MDF
mdf2 = ri.MDF(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.7174, 
              x0=103, A=20.0, f1=0.7, f2=-2, common_scale=1, 
              R_HM=2.2, X_h=20.4, psi=0.1571)
#mdf2.fit_edp()
#print("Report on MDF intensity")
#mdf2.report_edp()  

###############################################################################
# Work on S1G
s1g2 = ri.S1G(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.7174,
              x0=103, A=18.5, 
              rho_H1=10.77, Z_H1=20.86, sigma_H1=3.43,
              rho_M=9.23, sigma_M=1.67, psi=0.0873, common_scale=1)
s1g2.edp_par['rho_H1'].vary = False
s1g2.edp_par['sigma_H1'].vary = False
s1g2.edp_par['rho_M'].vary = False
s1g2.edp_par['sigma_M'].vary = False 
#s1g2.fit_edp()
#print("Report on S1G intensity")
#s1g2.report_edp()            

###############################################################################
# Work on M1G
m1g2 = ri.M1G(h, k, q, I, sigma, D=57.94, lambda_r=141.7, gamma=1.7174,
              x0=103, A=19.0, f1=0.6, f2=-1, 
              rho_H1=10.77, Z_H1=19.3, sigma_H1=3.43,
              rho_M=9.23, sigma_M=1.67, psi=0.157, common_scale=1)
m1g2.edp_par['rho_H1'].vary = False
m1g2.edp_par['sigma_H1'].vary = False
m1g2.edp_par['rho_M'].vary = False
m1g2.edp_par['sigma_M'].vary = False 
m1g2.fit_edp()
print("Report on M1G intensity")
m1g2.report_edp()   
