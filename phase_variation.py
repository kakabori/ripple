###############################################################################
def flatness(x):
    """Compute the flatness of a series of points"""
    chisqr = ((x-x.mean())**2).sum() / (x.size-1)
    return chisqr
    
def most_flat_profile(rip):
    mydict = {}
    best = 10**10
    # Grab indices for h = 6 orders
    #index = np.where((rip.h==6)|(rip.h==7)|(rip.h==9))
    index = np.where((rip.h==5))
    # There are N possible combinations for the phase factors
    N = 2**len(index[0])
    # Basic idea: convert i to binary string s. This loop will go through
    # all N possible combinations. 
    for i in xrange(N):
        s = generate_binary_string(i, N-1)
        a = np.array(binary_string_to_list(s), int)
        # Change binary 0 to phase factor -1
        a[a==0] = -1
        # Replace a subset of the phase factors with a generated combination
        rip.phase[index] = a       
        X, Z, DIST, ED = rip.get_EDP_between_two_points(start=(-40,21), 
                                                        end=(40,38), N=161)
        if flatness(ED) < best:
            best = flatness(ED)
            best_array = a
        mydict[s] = flatness(ED)
        #mylist.append(flatness(ED))
    rip.phase[index] = best_array
    return mydict, best_array, best
    
def binary_string_to_list(s):
    mylist = []
    for i in s:
        mylist.append(int(i))
    return mylist
    
def generate_binary_string(n, N):
    ret = "{0:b}".format(n)
    while len(ret) < len("{0:b}".format(N)):
        ret = '0' + ret
    return ret
  
def symmetry(x):
    symm = 0
    center = np.argmin(x)
    left_side = x[center-1::-1]
    right_side = x[center+1:]
    for i, j in zip(left_side, right_side):
        symm = symm + (i - j)**2
    return symm
    
def most_symmetric_profile(rip, val):
    mydict = {}
    best = 10**10
    # Grab indices 
    #index = np.where((rip.h==6)|(rip.h==7)|(rip.h==9))
    index = np.where((rip.h==val))
    # There are N possible combinations for the phase factors
    N = 2**len(index[0])
    # Basic idea: convert i to binary string s. This loop will go through
    # all N possible combinations. 
    for i in xrange(N):
        s = generate_binary_string(i, N-1)
        a = np.array(binary_string_to_list(s), int)
        # Change binary 0 to phase factor -1
        a[a==0] = -1
        # Replace a subset of the phase factors with a generated combination
        rip.phase[index] = a       
        X, Z, DIST, ED = rip.get_EDP_angle(center=(10,2.1), angle=-11.8, 
                                           length=100, stepsize=0.1)
        if symmetry(ED) < best:
            best = symmetry(ED)
            best_array = a
        mydict[s] = symmetry(ED)
    rip.phase[index] = best_array
    return mydict, best_array, best
    
###############################################################################
def linear(s=1, s_sigma=0.1, i=1, i_sigma=0.1):
    N = 1000
    x = np.linspace(-50,50,101)
    y = np.zeros((N,101))
    for index in range(N):
        slope = random.gauss(s, s_sigma)
        intercept = random.gauss(i, i_sigma)
        y[index] = slope*x + intercept
    mean = y.mean(axis=0)
    std = y.std(axis=0)
    plt.errorbar(x, mean, yerr=std, linestyle="None")

def resample_EDP(r, num=10000):
    """
    Monte-Carlo-resample the ripple EDP.
    
    Given form factors and corresponding uncertainties,
    this method randomly resamples the form factors averaging about
    the existing form factors, following Gaussian distributions
    with their sigma equal to the existing uncertainties.
    
    Input
    =====
    r: ripple object.
    N: number of resampling. The higher N is, the smoother the result is.
    """
    origin = (10, 2.1)
    angle = -11.8
    length = 100
    stepsize = 0.1
    x, z = origin
    N = length/stepsize + 1
    angle = angle*pi/180 
    if angle==0:
        # If angle is zero, the slope is infinite. 
        # In this case, x is constant.
        xpoints = x * np.ones(N)
        zpoints = np.linspace(z-length/2, z+length/2, N)
    else:
        slope = 1 / tan(angle)
        intercept = z - slope*x
        xpoints = np.linspace(x-length*sin(angle)/2, x+length*sin(angle)/2, N)
        zpoints = slope * xpoints + intercept  
        
    mylist = []      
    F_resample = generate_gauss_array(r.phase*r.F, r.sigma_F, num)
    for F in F_resample:
        DIST, EDP = calc_EDP(xpoints, zpoints, origin, r.qx, r.qz, F)
        mylist.append(EDP)
    mylist = np.array(mylist)
    mean = mylist.mean(axis=0)
    std = mylist.std(axis=0)
    plt.errorbar(DIST, mean, yerr=std, linestyle="None")        
    
def calc_EDP(xpoints, zpoints, origin, qx, qz, F):
    """Calculate 1D electron density profile from two dimensional 
    diffraction peak intensity data, at points specified by 
    xpoints and zpoints arrays.
    
    qx, qz, and F must be numpy arrays.
    """
    x0, z0 = origin
    EDP = []
    DIST = []
    for x, z in zip(xpoints, zpoints):
        tmp = F * np.cos(qx*x+qz*z)
        dist = np.sign(z-z0)*np.sqrt((x-x0)**2 + (z-z0)**2)
        DIST.append(dist)
        EDP.append(tmp.sum())  
    return DIST, EDP 

def generate_gauss_array(F, sigma, N=10000):
    """Generate an array whose element follows a normal distribution
    
    Return object is a 2D numpy arrays with each row 
    corresponding to an instance of resampled input F.
    
    example
    =======
    If input arrays are F = [1,20,300] and sigma=[0.1,2,30],
    returned list will look something like
    [[1.1, 20, 315],
     [0.9, 21, 299],
     [1.2, 18, 342],
     ...]
    """
    tmp = []
    for f, s in zip(F, sigma):
        tmp.append(np.random.normal(f, s, N))
    tmp = np.array(tmp)
    return np.transpose(tmp)
    
def generate_phase_varied_EDP(rip):
    mylist = []
    index = np.where((rip.h==6)|(rip.h==9))
    # There are N possible combinations for the phase factors
    N = 2**len(index[0])
    # Basic idea: convert i to binary string s. This loop will go through
    # all N possible combinations. 
    for i in xrange(N):
        s = generate_binary_string(i, N-1)
        a = np.array(binary_string_to_list(s), int)
        # Change binary 0 to phase factor -1
        a[a==0] = -1
        # Replace a subset of the phase factors with a generated combination
        rip.phase[index] = a       
        X, Z, DIST, ED = rip.get_EDP_between_two_points(start=(-40,21), 
                                                        end=(40,38), N=161)
        mylist.append(ED)
    mylist = np.array(mylist)
    mean = mylist.mean(axis=0)
    std = mylist.std(axis=0)
    return mean, std
