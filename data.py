import numpy as np

class Data(object):
    """This class implements a convenient access to X-ray form factor data 
    points by allowing a user to specify which data points are masked. 
    
    All variables are numpy objects with the same length.
    h, k : Miller indices
    qx, qz : qx and qz values
    F : signed form factor
    mask: specify which data points are returned in method calls
    
    For the mask variable, it accepts boolean. When data points 
    are set to True, returned arrays such as F exclude those points.
    """
    def __init__(self, h=None, k=None, qx=None, qz=None, F=None):
        """Inputs should be lists. They are stored as numpy arrays.
        
        h, k : Miller indices
        qx, qz : corresponding qx and qz values
        F : signed form factors
        """
        self.update(h, k, qx, qz, F)   
    
    def update(self, h=None, k=None, qx=None, qz=None, F=None):
        if h is not None:
            self.h = np.array(h, int)
            self.mask = np.zeros(len(self.h), bool)
        if k is not None:
            self.k = np.array(k, int)
        if qx is not None:
            self.qx = np.array(qx, float)
        if qz is not None:
            self.qz = np.array(qz, float)
        if F is not None:
            self.F = np.array(F, float) 
                        
    def form_factors(self, use_mask=True):
        if use_mask is True:
            F = self.F[~self.mask]
        else:
            F = self.F 
        return F
           
    def Miller_indices(self, use_mask=True):
        if use_mask is True:
            h = self.h[~self.mask]
            k = self.k[~self.mask]
        else:
            h = self.h
            k = self.k
        return h, k
    
    def qx_qz(self, use_mask=True):
        if use_mask is True:
            qx = self.qx[~self.mask]
            qz = self.qz[~self.mask] 
        else:
            qx = self.qx
            qz = self.qz
        return qx, qz        
               
    def set_mask(self, value, h, k=None):
        if k is None:
            self.mask[(self.h==h)] = value
        else:
            self.mask[(self.h==h)&(self.k==k)] = value
        
    def show_data(self):
        for h, k, F, mask in zip(self.h, self.k, self.F, self.mask):
            print(h, k, F, mask)
            
    def flip_phases(self, hin, kin):
        """Input can be single integer or lists"""
        for h, k in zip(hin, kin):
            self.F[(self.h==h)&(self.k==k)] *= -1
