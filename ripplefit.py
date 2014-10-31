from ripplefit.user_functions import *



if __name__ == '__main__':
    print("""
Typical commands:

load_data(filename="yourintensitydata.dat")

fit_edp()

For plot and export files, see user_manual.txt.
""")

    infilename = 'intensity/085_h9_ver6.dat'
    load_data(infilename)    
