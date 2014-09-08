#!/usr/bin/python

from Tkinter import *
import IPython 
from ripplefit import *

class StdoutRedirector(object):
    def __init__(self, text_widget):
        self.text_space = text_widget
        
    def write(self, string):
        self.text_space.insert('end', string)
        self.text_space.see('end')
        
class ripplefitGUI(object):
    def __init__(self, parent):
        infilename = 'intensity/085_h9_ver4.dat'
        h, k, q, I, sigma, combined = read_data_5_columns(infilename)
        self.m2g = M2G(h, k, q, I, sigma, D=57.8, lambda_r=145.0, gamma=1.714,
                       xM=90, A=25, f1=1.5, f2=-20, 
                       rho_H1=9.91, Z_H1=20, sigma_H1=2.94,
                       rho_H2=7.27, Z_H2=20, sigma_H2=1.47, 
                       rho_M=10.91, sigma_M=1.83, psi=0.1, common_scale=3)
        self.m2g.edp_par['f2'].vary = True
        self.m2g.edp_par['rho_H1'].vary = False
        self.m2g.edp_par['sigma_H1'].vary = False
        self.m2g.edp_par['rho_H2'].vary = False
        self.m2g.edp_par['sigma_H2'].vary = False
        self.m2g.edp_par['rho_M'].vary = False
        self.m2g.edp_par['sigma_M'].vary = False
        
        self.parent = parent
        self.initUI()
        
    def initUI(self):
        self.row_counter = {}
        self.row_counter['lattice'] = 0
        self.row_counter['contour'] = 0
        self.row_counter['transbilayer'] = 0
        
        top = Frame(self.parent)
        top.pack(side='top')
        bottom = Frame(self.parent)
        bottom.pack(side='top')
        left = Frame(top)
        left.pack(side='left')
        middle = Frame(top)
        middle.pack(side='left')
        right = Frame(top)
        right.pack(side='left')
        right2 = Frame(top)
        right2.pack(side='left')

        self.lattice_pars = {}    
        self.lattice_pars['D'] = DoubleVar()
        self.lattice_pars['D'].set(57.8)
        self.lattice_pars['lambda_r'] = DoubleVar()
        self.lattice_pars['lambda_r'].set(145.0)
        self.lattice_pars['gamma'] = DoubleVar()
        self.lattice_pars['gamma'].set(1.714)
        self.contour_pars = {}
        self.contour_pars['A'] = DoubleVar()
        self.contour_pars['A'].set(22)
        self.contour_pars['xM'] = DoubleVar()
        self.contour_pars['xM'].set(90)
        self.contour_pars['f1'] = DoubleVar()
        self.contour_pars['f1'].set(1.5)
        self.contour_pars['f2'] = DoubleVar()
        self.contour_pars['f2'].set(-20)
        self.transbilayer_pars = {}
        self.transbilayer_pars['rho_H1'] = DoubleVar()
        self.transbilayer_pars['rho_H1'].set(9.91)
        self.transbilayer_pars['Z_H1'] = DoubleVar()
        self.transbilayer_pars['Z_H1'].set(20)
        self.transbilayer_pars['sigma_H1'] = DoubleVar()
        self.transbilayer_pars['sigma_H1'].set(2.94)
        self.transbilayer_pars['rho_H2'] = DoubleVar()
        self.transbilayer_pars['rho_H2'].set(7.27)
        self.transbilayer_pars['Z_H2'] = DoubleVar()
        self.transbilayer_pars['Z_H2'].set(20)
        self.transbilayer_pars['sigma_H2'] = DoubleVar()
        self.transbilayer_pars['sigma_H2'].set(1.47)
        self.transbilayer_pars['rho_M'] = DoubleVar()
        self.transbilayer_pars['rho_M'].set(10.91)
        self.transbilayer_pars['sigma_M'] = DoubleVar()
        self.transbilayer_pars['sigma_M'].set(1.83)
        self.transbilayer_pars['psi'] = DoubleVar()
        self.transbilayer_pars['psi'].set(0.1)
        self.transbilayer_pars['common_scale'] = DoubleVar()
        self.transbilayer_pars['common_scale'].set(3)

        for key in ['D', 'lambda_r', 'gamma']:
            #lattice_pars[key].trace("w", update_lattice_pars)
            tmp, self.row_counter['lattice'] = \
                self.textentry(left, self.lattice_pars[key], key, 
                               self.row_counter['lattice'])

        for key in ['A', 'xM', 'f1', 'f2']:
            #contour_pars[key].trace("w", update_contour_pars)
            tmp, self.row_counter['contour'] = \
                self.textentry(middle, self.contour_pars[key], key,
                               self.row_counter['contour'])

        for key in ['rho_H1', 'Z_H1', 'sigma_H1', 'rho_H2', 'Z_H2', 'sigma_H2', 
                    'rho_M', 'sigma_M', 'psi', 'common_scale']:
            #transbilayer_pars[key].trace("w", update_transbilayer_pars)    
            tmp, self.row_counter['transbilayer'] = \
                self.textentry(right, self.transbilayer_pars[key], key,
                               self.row_counter['transbilayer'])
              
        #self.text_box = Text(right2, wrap=None, height=11, width=50)
        #self.text_box.grid(column=0, row=0, columnspan=2, sticky='NSWE', 
        #                   padx=5, pady=5)
        #sys.stdout = StdoutRedirector(self.text_box)
        
        fit_btn = Button(bottom, text='Fit', command=self.fit_model)
        fit_btn.pack(side='left')
        
        # Make a quit button
        quit_button = Button(bottom, text='quit', command=self.quit,
                             background='yellow', foreground='blue')
        quit_button.pack(side='left', pady=5, fill='x')
        #self.master.bind('<q>', self.quit)
        
    def textentry(self, parent, variable, label, row_counter):
        L = Label(parent, text=label)
        L.grid(column=0, row=row_counter, sticky='w')
        widget = Entry(parent, textvariable=variable, width=8)
        widget.grid(column=1, row=row_counter)
        row_counter += 1
        return widget, row_counter
    
    def quit(self, event=None):
        self.parent.destroy()
    
    def fit_model(self):
        for key in self.lattice_pars:
            self.m2g.latt_par[key].value = self.lattice_pars[key].get()
        for key in self.contour_pars:
            self.m2g.edp_par[key].value = self.contour_pars[key].get()
        for key in self.transbilayer_pars:
            self.m2g.edp_par[key].value = self.transbilayer_pars[key].get()
        self.m2g.fit_edp()
        self.m2g.report_edp()
        for key in self.contour_pars:
            self.contour_pars[key].set(self.m2g.edp_par[key].value)
        for key in self.transbilayer_pars:
            self.transbilayer_pars[key].set(self.m2g.edp_par[key].value)


    
    
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


       





root = Tk()
ripple = ripplefitGUI(root)
IPython.embed()
#root.mainloop()

