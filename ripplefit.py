import ripple as r
from user_functions import *

print("""
Typical commands:

load_data(filename="yourintensitydata.dat")

edp_fit()

For plot and export files, see user_manual.txt.
""")

infilename = 'intensity/085_h9_ver6.dat'
load_data(infilename)
