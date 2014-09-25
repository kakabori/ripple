#!/usr/bin/python

from Tkinter import *
#from ttk import *
#import ttk
import IPython 
from ripplefit import *
import Pmw
import tkMessageBox, tkFileDialog, tkColorChooser
import string, sys, os

class StdoutRedirector(object):
    def __init__(self, text_widget):
        self.text_space = text_widget
        
    def write(self, string):
        self.text_space.insert('end', string)
        self.text_space.see('end')
        
class ripplefitGUI(object):
    def __init__(self, parent):
        infilename = 'intensity/085_h9_ver6.dat'
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
        
        self.h = h
        self.k = k
        self.I = I
        self.sigma = sigma
        
        self.parent = parent
        self.pulldown_menus(self.parent)
        top = Frame(self.parent)
        top.pack(side='top')
        self.pf = ParameterFields(top)
        
        bottom = Frame(self.parent)
        bottom.pack(side='top')
        fit_btn = Button(bottom, text='Fit', command=self.fit_model)
        fit_btn.pack(side='left')
               
    def pulldown_menus(self, parent):
        self.menu_bar = Pmw.MenuBar(parent,
                                    hull_relief='raised',
                                    hull_borderwidth=1)
        self.menu_bar.pack(fill='x')

        self.menu_bar.addmenu('File', None, tearoff=True)
        self.menu_bar.addmenuitem('File', 'command', 
             statusHelp='Open a file',
             label='Open...', 
             command=self.file_read)
        self.menu_bar.addmenuitem('File', 'command', 
             statusHelp='Save a file',
             label='Save as...', 
             command=self.file_save)
        self.menu_bar.addmenuitem('File', 'command', 
             statusHelp='Exit this application',
             label='Quit',
             command=self.quit)

        self.menu_bar.addmenu('Dialogs', None, tearoff=True)
        self.menu_bar.addmenuitem('Dialogs', 'command',
             label='Pmw user-defined dialog',
             command=self.userdef_dialog)

    def file_read(self):
        fname = tkFileDialog.Open(filetypes=[('anyfile','*')]).show()
        text = 'chosen file to open: ' + os.path.basename(fname)
        # the dialog checks the validity of the filename, but
        # pressing Cancel results in an empty return string
        if fname:
            self.display_file(fname, self.master)

    def file_save(self):
        fname = tkFileDialog.SaveAs(
                filetypes=[('temporary files','*.tmp')],
                initialfile='myfile.tmp',
                title='Save a file').show()
        text = 'chosen file to save: "' + os.path.basename(fname) + '"'
        
    def quit(self, event=None):
        self.parent.destroy()
        
    def userdef_dialog(self):
        self.userdef_d = Pmw.Dialog(self.parent,
                          title='Programmer-Defined Dialog',
                          buttons=('Update', 'Cancel'),
                          #defaultbutton='Apply',
                          command=self.userdef_dialog_action)

        self.userdef_d_gui = InputDataFields(self.userdef_d.interior(),
                                             self.h, self.k, self.I, self.sigma)
        self.userdef_d_gui.pack()

    def userdef_dialog_action(self, result):
        # result contains the name of the button that we clicked
        if result == 'Apply':
            # example on extracting dialog variables:
            case = self.userdef_d_gui.case.get()
        else:
            text = 'you just canceled the dialog'
        self.userdef_d.destroy()  # destroy dialog window   
        
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
            
                    
class ParameterFields(object):
    """
    Create a GUI for parameters
    """
    def __init__(self, parent):
        self.parent = parent
        left_frame = Frame(self.parent)
        #left_frame.pack(side='left')
        left_frame.grid(column=0, row=0, sticky='n')
        right_frame = Frame(self.parent)
        #right_frame.pack(side='right')
        right_frame.grid(column=1, row=0)
        self.row_counter = {}
        self.row_counter['left'] = 0
        self.row_counter['right'] = 0

        # Initialize associated variables
        # Call as lattice_pars['D']['val'].get() etc.
        self.lattice_pars = {}
        for key, i in zip(['D', 'lambda_r', 'gamma'], [57.8, 145.0, 1.714]):
            self.lattice_pars[key] = {'val': DoubleVar(), 'vary': IntVar()}    
            self.lattice_pars[key]['val'].set(i)
            self.lattice_pars[key]['vary'].set(False)
        for key in ['D', 'lambda_r', 'gamma']:
            tmp, self.row_counter['left'] = \
                self.textentry(left_frame, self.lattice_pars[key]['val'], 
                               self.lattice_pars[key]['vary'], key, 
                               self.row_counter['left'])
                                           
        self.contour_pars = {}
        for key, i in zip(['A', 'xM', 'f1', 'f2'], [22.0, 90.0, 1.5, -20.0]):
            self.contour_pars[key] = {'val': DoubleVar(), 'vary': IntVar()}
            self.contour_pars[key]['val'].set(i)
            self.contour_pars[key]['vary'].set(True)
        for key in ['A', 'xM', 'f1', 'f2']:
            tmp, self.row_counter['left'] = \
                self.textentry(left_frame, self.contour_pars[key]['val'], 
                               self.contour_pars[key]['vary'], key,
                               self.row_counter['left'])
            
        self.transbilayer_pars = {}
        for key, i in zip(['rho_H1', 'Z_H1', 'sigma_H1', 'rho_H2', 'Z_H2', 
                           'sigma_H2', 'rho_M', 'sigma_M', 'psi', 
                           'common_scale'], [9.94, 20.0, 2.94, 7.27, 20.0, 
                           1.47, 10.91, 1.83, 0.1, 3]):
            self.transbilayer_pars[key] = {'val': DoubleVar(), 'vary': IntVar()}
            self.transbilayer_pars[key]['val'].set(i)
            self.transbilayer_pars[key]['vary'].set(True)
        for key in ['rho_H1', 'Z_H1', 'sigma_H1', 'rho_H2', 'Z_H2', 'sigma_H2', 
                    'rho_M', 'sigma_M', 'psi', 'common_scale']:  
            tmp, self.row_counter['right'] = \
                self.textentry(right_frame, self.transbilayer_pars[key]['val'], 
                               self.transbilayer_pars[key]['vary'], key,
                               self.row_counter['right'])       
        
    def textentry(self, parent, text_var, check_var, label, row_counter):
        C = Checkbutton(parent, text=label, variable=check_var)
        C.grid(column=0, row=row_counter, sticky='w')
        widget = Entry(parent, textvariable=text_var, width=8)
        widget.grid(column=1, row=row_counter)
        row_counter += 1
        return widget, row_counter

               
class InputDataFields(object):
    """
    Create a GUI for input data points
    """
    def __init__(self, parent, h, k, I, sigma): 
        """
        Create widgets.
        
        parent    parent widget
        h         h list
        k         k list
        I         Intensity list
        sigma     uncertainty list
        """
        self.parent = parent
        #self.topframe = Pmw.ScrolledFrame(self.parent, usehullsize=1, 
        #                                  hull_height=210, hull_width=340)
        self.topframe = Pmw.ScrolledFrame(self.parent)
        self.create(self.topframe.interior(), h, k, I, sigma)

    def pack(self, **kwargs):
        """
        Pack the topframe. The location of the InputFields GUI in
        the parent widget can be controlled by the calling code.
        """
        self.topframe.pack(kwargs, expand=True, fill='both')
        
    def create(self, parent, hL, kL, IL, sigmaL):
        """Create all widgets."""
        row_counter = 0
        h_widgets, k_widgets, I_widgets, sigma_widgets = [], [], [], []
        for h, k, I, sigma in zip(hL, kL, IL, sigmaL):
            widget = Entry(parent, width=8)
            widget.insert('end', str(h))
            widget.grid(column=0, row=row_counter)
            h_widgets.append(widget)
            widget = Entry(parent, width=8)
            widget.insert('end', str(k))
            widget.grid(column=1, row=row_counter)
            k_widgets.append(widget)
            widget = Entry(parent, width=8)
            widget.insert('end', str(I))
            widget.grid(column=2, row=row_counter)      
            I_widgets.append(widget)
            widget = Entry(parent, width=8)
            widget.insert('end', str(sigma))
            widget.grid(column=3, row=row_counter)
            sigma_widgets.append(widget)
            row_counter = row_counter + 1      

  

if __name__ == '__main__':
    root = Tk()
    #h=[1,1,1]
    #k=[-1,0,1]
    #I=[73,102,34]
    #sigma=[0.1,3,4]
    #a = InputDataFields(root, h, k, I, sigma)
    #a.pack()
    ripple = ripplefitGUI(root)
    #IPython.embed()
    root.mainloop()

