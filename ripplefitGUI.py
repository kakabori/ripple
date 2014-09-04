#!/usr/bin/python

from Tkinter import *
import IPython 
from ripplefit import *

infilename = 'intensity/085_h9_ver4.dat'
h, k, q, I, sigma, combined = read_data_5_columns(infilename)
m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
          xM=90, A=25, f1=1.5, f2=-20, 
          rho_H1=9.91, Z_H1=20, sigma_H1=2.94,
          rho_H2=7.27, Z_H2=20, sigma_H2=1.47, 
          rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=3)
m2g.edp_par['f2'].vary = True
m2g.edp_par['rho_H1'].vary = False
m2g.edp_par['sigma_H1'].vary = False
m2g.edp_par['rho_H2'].vary = False
m2g.edp_par['sigma_H2'].vary = False
m2g.edp_par['rho_M'].vary = False
m2g.edp_par['sigma_M'].vary = False

root = Tk()
top = Frame(root)
top.pack(side='top')
bottom = Frame(root)
bottom.pack(side='top')
left = Frame(top)
left.pack(side='left')
middle = Frame(top)
middle.pack(side='left')
right = Frame(top)
right.pack(side='left')

def textentry(parent, variable, label):
    f = Frame(parent)
    f.pack(side='top', padx=2, pady=2)
    L = Label(f, text=label)
    L.pack(side='left')
    widget = Entry(f, textvariable=variable, width=8)
    widget.pack(side='left', anchor='w')
    return widget

lattice_pars = {}    
lattice_pars['D'] = DoubleVar()
lattice_pars['D'].set(57.8)
lattice_pars['lambda_r'] = DoubleVar()
lattice_pars['lambda_r'].set(145.0)
lattice_pars['gamma'] = DoubleVar()
lattice_pars['gamma'].set(1.714)
contour_pars = {}
contour_pars['A'] = DoubleVar()
contour_pars['A'].set(22)
contour_pars['xM'] = DoubleVar()
contour_pars['xM'].set(90)
contour_pars['f1'] = DoubleVar()
contour_pars['f1'].set(1.5)
contour_pars['f2'] = DoubleVar()
contour_pars['f2'].set(-20)
transbilayer_pars = {}
transbilayer_pars['rho_H1'] = DoubleVar()
transbilayer_pars['rho_H1'].set(9.91)
transbilayer_pars['Z_H1'] = DoubleVar()
transbilayer_pars['Z_H1'].set(20)
transbilayer_pars['sigma_H1'] = DoubleVar()
transbilayer_pars['sigma_H1'].set(2.94)
transbilayer_pars['rho_H2'] = DoubleVar()
transbilayer_pars['rho_H2'].set(7.27)
transbilayer_pars['Z_H2'] = DoubleVar()
transbilayer_pars['Z_H2'].set(20)
transbilayer_pars['sigma_H2'] = DoubleVar()
transbilayer_pars['sigma_H2'].set(1.47)
transbilayer_pars['rho_M'] = DoubleVar()
transbilayer_pars['rho_M'].set(10.91)
transbilayer_pars['sigma_M'] = DoubleVar()
transbilayer_pars['sigma_M'].set(1.83)
transbilayer_pars['psi'] = DoubleVar()
transbilayer_pars['psi'].set(0.1)
transbilayer_pars['common_scale'] = DoubleVar()
transbilayer_pars['common_scale'].set(3)

def update_lattice_pars(*args):
    global lattice_pars
    for key in lattice_pars:
        m2g.latt_par[key].value = lattice_pars[key].get()

def update_contour_pars(*args):
    global contour_pars
    for key in contour_pars:
        m2g.edp_par[key].value = contour_pars[key].get()
        
def update_transbilayer_pars(*args):
    global transbilayer_pars
    for key in transbilayer_pars:
        m2g.edp_par[key].value = transbilayer_pars[key].get()        

for key in lattice_pars:
    #lattice_pars[key].trace("w", update_lattice_pars)
    textentry(left, lattice_pars[key], key)
    
for key in contour_pars:
    #contour_pars[key].trace("w", update_contour_pars)
    textentry(middle, contour_pars[key], key)
    
for key in transbilayer_pars:
    #transbilayer_pars[key].trace("w", update_transbilayer_pars)    
    textentry(right, transbilayer_pars[key], key)

def fit_model():
    for key in lattice_pars:
        m2g.latt_par[key].value = lattice_pars[key].get()
    for key in contour_pars:
        m2g.edp_par[key].value = contour_pars[key].get()
    for key in transbilayer_pars:
        m2g.edp_par[key].value = transbilayer_pars[key].get()
    m2g.fit_edp()
    m2g.report_edp()
    for key in contour_pars:
        contour_pars[key].set(m2g.edp_par[key].value)
    for key in transbilayer_pars:
        transbilayer_pars[key].set(m2g.edp_par[key].value)

fit_btn = Button(bottom, text='Fit', command=fit_model)
fit_btn.pack(side='top')

IPython.embed()
#root.mainloop()

