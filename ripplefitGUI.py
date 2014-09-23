#!/usr/bin/python

from Tkinter import *
#from ttk import *
import ttk
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

        # Initialize associated variables
        # Call as lattice_pars['D']['val'].get() etc.
        self.lattice_pars = {}
        for key, i in zip(['D', 'lambda_r', 'gamma'], [57.8, 145.0, 1.714]):
            self.lattice_pars[key] = {'val': DoubleVar(), 'vary': IntVar()}    
            self.lattice_pars[key]['val'].set(i)
            self.lattice_pars[key]['vary'].set(False)
            
        self.contour_pars = {}
        for key, i in zip(['A', 'xM', 'f1', 'f2'], [22.0, 90.0, 1.5, -20.0]):
            self.contour_pars[key] = {'val': DoubleVar(), 'vary': IntVar()}
            self.contour_pars[key]['val'].set(i)
            self.contour_pars[key]['vary'].set(True)
            
        self.transbilayer_pars = {}
        for key, i in zip(['rho_H1', 'Z_H1', 'sigma_H1', 'rho_H2', 'Z_H2', 
                           'sigma_H2', 'rho_M', 'sigma_M', 'psi', 
                           'common_scale'], [9.94, 20.0, 2.94, 7.27, 20.0, 
                           1.47, 10.91, 1.83, 0.1, 3]):
            self.transbilayer_pars[key] = {'val': DoubleVar(), 'vary': IntVar()}
            self.transbilayer_pars[key]['val'].set(i)
            self.transbilayer_pars[key]['vary'].set(True)
                    
        for key in ['D', 'lambda_r', 'gamma']:
            #lattice_pars[key].trace("w", update_lattice_pars)
            tmp, self.row_counter['lattice'] = \
                self.textentry(left, self.lattice_pars[key]['val'], 
                               self.lattice_pars[key]['vary'], key, 
                               self.row_counter['lattice'])

        for key in ['A', 'xM', 'f1', 'f2']:
            #contour_pars[key].trace("w", update_contour_pars)
            tmp, self.row_counter['contour'] = \
                self.textentry(middle, self.contour_pars[key]['val'], 
                               self.contour_pars[key]['vary'], key,
                               self.row_counter['contour'])

        for key in ['rho_H1', 'Z_H1', 'sigma_H1', 'rho_H2', 'Z_H2', 'sigma_H2', 
                    'rho_M', 'sigma_M', 'psi', 'common_scale']:
            #transbilayer_pars[key].trace("w", update_transbilayer_pars)    
            tmp, self.row_counter['transbilayer'] = \
                self.textentry(right, self.transbilayer_pars[key]['val'], 
                               self.transbilayer_pars[key]['vary'], key,
                               self.row_counter['transbilayer'])
        
        fit_btn = Button(bottom, text='Fit', command=self.fit_model)
        fit_btn.pack(side='left')
        
        # Make a quit button
        quit_button = Button(bottom, text='quit', command=self.quit,
                             background='yellow', foreground='blue')
        quit_button.pack(side='left', pady=5, fill='x')
        #self.master.bind('<q>', self.quit)
        
        n = ttk.Notebook(bottom)
        f1 = ttk.Frame(n); # first page, which would get widgets gridded into it
        f2 = ttk.Frame(n); # second page
        n.add(f1, text='One')
        n.add(f2, text='Two')
        
    def textentry(self, parent, text_var, check_var, label, row_counter):
        #L = Label(parent, text=label)
        #L.grid(column=0, row=row_counter, sticky='w')
        C = Checkbutton(parent, text=label, variable=check_var)
        C.grid(column=0, row=row_counter, sticky='w')
        widget = Entry(parent, textvariable=text_var, width=8)
        widget.grid(column=1, row=row_counter)
        row_counter += 1
        return widget, row_counter
    
    def quit(self, event=None):
        self.parent.destroy()
    
    def fit_model(self):
        """
        Update the M2G parameter objects with the current values in the GUI,
        perform a NLSQ optimization, and then update the GUI with the best
        fit values.
        """
        # Update Parameter objects with values and free/fix from the GUI
        for key in self.lattice_pars:
            self.m2g.latt_par[key].value = self.lattice_pars[key]['val'].get()
        for key in self.contour_pars:
            self.m2g.edp_par[key].value = self.contour_pars[key]['val'].get()
            self.m2g.edp_par[key].vary = bool(self.contour_pars[key]['vary'].get())
        for key in self.transbilayer_pars:
            self.m2g.edp_par[key].value = self.transbilayer_pars[key]['val'].get()
            self.m2g.edp_par[key].vary = bool(self.transbilayer_pars[key]['vary'].get())                     
        self.m2g.fit_edp()
        self.m2g.report_edp()
        for key in self.contour_pars:
            self.contour_pars[key]['val'].set(self.m2g.edp_par[key].value)
        for key in self.transbilayer_pars:
            self.transbilayer_pars[key]['val'].set(self.m2g.edp_par[key].value)





       





root = Tk()
ripple = ripplefitGUI(root)
IPython.embed()
#root.mainloop()

