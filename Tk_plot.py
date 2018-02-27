#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 11:34:32 2018

@author: ryanday
"""

import matplotlib
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'figure.autolayout':True})
matplotlib.use('TkAgg')

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from matplotlib.figure import Figure
import matplotlib.cm as cm

import sys


if sys.version_info[0]<3:
    import Tkinter as Tk
else:
    import tkinter as Tk
    


class plot_intensity_interface:
    
    def __init__(self,expmnt,Adict):
        self.root = Tk.Tk()
        self.expmnt = expmnt
        self.Adict = Adict
        print('Initializing spectral function...')
        _,self.Imat = self.expmnt.spectral(self.Adict)
        self.Imat_dict = {'I_0':self.Imat}
        self.x = self.Adict['cube']['X']
        self.y = self.Adict['cube']['Y']
        self.w = self.Adict['cube']['E']
        self.dx,self.dy,self.dw = (self.x[1]-self.x[0])/self.x[2],(self.y[1]-self.y[0])/self.y[2],(self.w[1]-self.w[0])/self.w[2]
        print('Initializing interface...')
        self.plot_make()

    def plot_make(self):
           
        self.root.wm_title('UBC TB-ARPES DATA PLOTTER')

        fig1 = Figure(figsize=(5,5))
        fig2 = Figure(figsize=(5,5))
        fig3 = Figure(figsize=(5,5))
        
        rc('font',**{'family':'serif','serif':['Palatino'],'size':12})
        rc('text',usetex = False) 

        ax1= fig1.add_subplot(111)
        ax2= fig2.add_subplot(111)
        ax3= fig3.add_subplot(111)
        sp1 = ax1.imshow(self.Imat[:,:,int(self.w[2]/2)],extent=[self.y[0],self.y[1],self.x[0],self.x[1]],cmap=cm.magma)
        sp2 = ax2.imshow(self.Imat[int(self.x[2]/2),:,:],extent=[self.w[0],self.w[1],self.y[0],self.y[1]],cmap=cm.magma)
        sp3 = ax3.imshow(self.Imat[:,int(self.y[2]/2),:],extent=[self.w[0],self.w[1],self.x[0],self.x[1]],cmap=cm.magma)
        sp1.set_clim(vmin=0,vmax=self.Imat.max()*1.05)
        sp2.set_clim(vmin=0,vmax=self.Imat.max()*1.05)
        sp3.set_clim(vmin=0,vmax=self.Imat.max()*1.05)
        asp1 = (self.x[1]-self.x[0])/(self.y[1]-self.y[0])
        asp2 = (self.w[1]-self.w[0])/(self.y[1]-self.y[0])
        asp3 = (self.w[1]-self.w[0])/(self.x[1]-self.x[0])
        ax1.set_aspect(asp1)
        ax2.set_aspect(asp2)
        ax3.set_aspect(asp3)
        
        ax1.tick_params(axis='both',labelsize=12)
        ax2.tick_params(axis='both',labelsize=12)
        ax3.tick_params(axis='both',labelsize=12)
        
        ax1.set_ylabel('Momentum y ($\AA^{-1}$)')
        ax1.set_xlabel('Momentum x ($\AA^{-1}$)')
        ax2.set_ylabel('Momentum y ($\AA^{-1}$)')
        ax2.set_xlabel('Energy (eV)')
        ax3.set_ylabel('Momentum x ($\AA^{-1}$)')
        ax3.set_xlabel('Energy (eV)')

        xpol_label = Tk.Label(master=self.root,text="Polarization x").grid(row=4,column=0)
        ypol_label = Tk.Label(master=self.root,text="Polarization y").grid(row=4,column=3)
        zpol_label = Tk.Label(master=self.root,text="Polarization z").grid(row=4,column=6)

        xpol_ent = Tk.Entry(master=self.root)
        ypol_ent = Tk.Entry(master=self.root)
        zpol_ent = Tk.Entry(master=self.root)

        xpol_ent.grid(row=4,column=1,columnspan=1)
        ypol_ent.grid(row=4,column=4,columnspan=1)
        zpol_ent.grid(row=4,column=7,columnspan=1)

        plt_win_1 = FigureCanvasTkAgg(fig1,master=self.root)
        plt_win_1.show()
        plt_win_1.get_tk_widget().grid(row=0,column=0,columnspan=3,rowspan=3)
        
        plt_win_2 = FigureCanvasTkAgg(fig2,master=self.root)
        plt_win_2.show()
        plt_win_2.get_tk_widget().grid(row=0,column=3,columnspan=3,rowspan=3)

        plt_win_3 = FigureCanvasTkAgg(fig3,master=self.root)
        plt_win_3.show()
        plt_win_3.get_tk_widget().grid(row=0,column=6,columnspan=3,rowspan=3)


        def _quit_win():
            self.root.quit()
            self.root.destroy()
    
        def _update_slice(event):
            sv1 = int((float(slide_1.get())-self.w[0])/(self.dw))
            sv2 = int((float(slide_2.get())-self.x[0])/(self.dx))
            sv3 = int((float(slide_3.get())-self.y[0])/(self.dy))
    
            sp1.set_data(self.Imat[:,:,sv1])
            sp2.set_data(self.Imat[sv2,:,:])
            sp3.set_data(self.Imat[:,sv3,:])
            fig1.canvas.draw()
            fig2.canvas.draw()
            fig3.canvas.draw()
    
        s1_label = Tk.Label(master=self.root,text="Energy (eV)").grid(row=3,column=0)
        slide_1 = Tk.Scale(master=self.root,from_=self.w[0],to=(self.w[1]-self.dw),orient='horizontal',resolution=self.dw,command=_update_slice)
        slide_1.grid(row=3,column=1,sticky='W')

        s2_label = Tk.Label(master=self.root,text="Momentum Kx (1/Å)").grid(row=3,column=3)
        slide_2 = Tk.Scale(master=self.root,from_=self.x[0],to=(self.x[1]-self.dx),orient='horizontal',resolution=self.dx,command=_update_slice)
        slide_2.grid(row=3,column=4,sticky='W')

        s3_label = Tk.Label(master=self.root,text="Momentum K_y (1/Å)").grid(row=3,column=6)
        slide_3 = Tk.Scale(master=self.root,from_=self.y[0],to=(self.y[1]-self.dy),orient='horizontal',resolution=self.dy,command=_update_slice)
        slide_3.grid(row=3,column=7,sticky='W')


        def _update_pol():

            try:
    
                self.Adict['pol'] = xyz(xpol_ent.get(),ypol_ent.get(),zpol_ent.get())
                print('Recalculating matrix elements...')
                _,self.Imat = self.expmnt.spectral(self.Adict)
                print('Matrix element calculation complete')
                
                sv1 = int((float(slide_1.get())-self.w[0])/(self.dw))
                sv2 = int((float(slide_2.get())-self.x[0])/(self.dx))
                sv3 = int((float(slide_3.get())-self.y[0])/(self.dy))
    
                sp1.set_data(self.Imat[:,:,sv1])
                sp2.set_data(self.Imat[sv2,:,:])
                sp3.set_data(self.Imat[:,sv3,:])
                sp1.set_clim(vmin=self.Imat.min(),vmax=self.Imat.max())
                sp2.set_clim(vmin=self.Imat.min(),vmax=self.Imat.max())
                sp3.set_clim(vmin=self.Imat.min(),vmax=self.Imat.max())
                fig1.canvas.draw()
                fig2.canvas.draw()
                fig3.canvas.draw()
            
            except ValueError:
                print('ERROR!! Polarization vector must be complex float!')
    
    
        update_pol_button = Tk.Button(master=self.root,text='Update Polarization',command=_update_pol)
        update_pol_button.grid(row=5,column=0)
        
        
        def _customize():
            
            cwin = Tk.Toplevel()
            cwin.wm_title('CUSTOMIZE MAP DISPLAY')
            #List of available Maps
            list_lab = Tk.Label(master=cwin,text='Available maps: ').grid(row=0,column=0)
            mat_list = (d for d in self.Imat_dict)
            mat_listbox = Tk.Listbox(master=cwin)
            mat_listbox.grid(row=1,column=0,columnspan=4)
            for d in list(enumerate(self.Imat_dict)):
                mat_listbox.insert("end",d[1])
            #Spin Vector
            spin_label=Tk.Label(master=cwin,text='Spin Axis: ').grid(row=2,column=0)
            sx = Tk.Entry(master=cwin)
            sy = Tk.Entry(master=cwin)
            sz = Tk.Entry(master=cwin)
            
            sx.grid(row=2,column=1)
            sy.grid(row=2,column=2)
            sz.grid(row=2,column=3)
            
            #Spin projection -- -1, None, 1
            proj_label = Tk.Label(master=cwin,text='Spin Projection: ').grid(row=3,column=0)
            proj_list = ("Down","Up","None")
            proj_var = Tk.StringVar()
            proj_var.set(proj_list[2])
            proj_opt = Tk.OptionMenu(cwin,proj_var,*proj_list)
            proj_opt.grid(row=3,column=1)
            #Polarization Vector
            pol_label = Tk.Label(master=cwin,text='Polarization: ').grid(row=4,column=0)
            pol_x = Tk.Entry(master=cwin)
            pol_y = Tk.Entry(master=cwin)
            pol_z = Tk.Entry(master=cwin)
            
            pol_x.grid(row=4,column=1)
            pol_y.grid(row=4,column=2)
            pol_z.grid(row=4,column=3)
            #Name of New Map -- option of string or Blank
            mat_nm_label = Tk.Label(master=cwin,text='Map Name: ').grid(row=5,column=0)
            mat_nm = Tk.Entry(master=cwin)
            mat_nm.grid(row=5,column=1)
        
            
            #Add a new map to the list of available maps
            def _add_map():
                
                tmp_dict = self.Adict.copy()
                #Define the polarization vector
                tmp_dict['pol'] = xyz(pol_x.get(),pol_y.get(),pol_z.get())
                
                proj_choice = proj_var.get()
                if proj_choice=="None":
                    tmp_dict['spin']=None
                else:
                    
                    tmp_dict['spin'] = [-1 if proj_choice=="Down" else 1,np.array([float(sx.get()),float(sy.get()),float(sz.get())])]
                print(tmp_dict['spin'])
                mat_name = mat_nm.get() if mat_nm.get()!="" else "I_{:d}".format(len(self.Imat_dict))

                _,self.Imat_dict[mat_name] = self.expmnt.spectral(tmp_dict)
                
                mat_listbox.insert("end",mat_name)            
            
            add_button = Tk.Button(master=cwin,text ='Add Map ',command=_add_map)
            add_button.grid(row=6,column=1,columnspan=2)            
            
            #add option of operating on datasets
            def _gen_map():
                st_raw = op_entry.get()
                for d in self.Imat_dict:
                    
                    replacement = 'self.Imat_dict["{:s}"]'.format(d)
                    st_raw = st_raw.replace(d,replacement)
                    
                    self.Imat_dict[d]+= abs(self.Imat_dict[d][np.nonzero(self.Imat_dict[d])]).min()*10**-4 #avoid divergence for division if zeros present
                st_raw = st_raw.replace("SQRT","np.sqrt")
                st_raw = st_raw.replace("COS","np.cos")
                st_raw = st_raw.replace("SIN","np.sin")
                st_raw = st_raw.replace("TAN","np.tan")
                st_raw = st_raw.replace("EXP","np.exp")
                st_raw = st_raw.replace("LOG","np.log")
                
                tmp_mat = eval(st_raw)
                map_nm = mat_nm.get() if (mat_nm.get()!="" or bool(sum([mat_nm==d for d in self.Imat_dict]))) else "I_{:d}".format(len(self.Imat_dict))
                mat_listbox.insert("end",map_nm)
                self.Imat_dict[map_nm] = tmp_mat
                
        
            op_label = Tk.Label(master=cwin,text="Operate: ").grid(row=7,column=0)
            op_entry = Tk.Entry(master=cwin)
            op_entry.grid(row=7,column=1)
            op_gen = Tk.Button(master=cwin,text = "Generate",command=_gen_map)
            op_gen.grid(row=7,column=2)
            
            def _plot_map():
                selected = mat_listbox.curselection()[0]
                try:
                    
                    self.Imat = self.Imat_dict[mat_listbox.get(selected)]
                    sv1 = int((float(slide_1.get())-self.w[0])/(self.dw))
                    sv2 = int((float(slide_2.get())-self.x[0])/(self.dx))
                    sv3 = int((float(slide_3.get())-self.y[0])/(self.dy))
    
                    sp1.set_data(self.Imat[:,:,sv1])
                    sp2.set_data(self.Imat[sv2,:,:])
                    sp3.set_data(self.Imat[:,sv3,:])
                    sp1.set_clim(vmin=self.Imat.min(),vmax=self.Imat.max())
                    sp2.set_clim(vmin=self.Imat.min(),vmax=self.Imat.max())
                    sp3.set_clim(vmin=self.Imat.min(),vmax=self.Imat.max())
                    fig1.canvas.draw()
                    fig2.canvas.draw()
                    fig3.canvas.draw()
                    cwin.quit()
                    cwin.destroy()
                    
                except IndexError:
                    print('Please select a map to plot before attempting to plot.')
            
            plot_button = Tk.Button(master=cwin,text="Plot Selection", command=_plot_map)
            plot_button.grid(row=8,column=0,columnspan=2)
            
            def _quit_sub():
                cwin.quit()
                cwin.destroy()
            
            close_button = Tk.Button(master= cwin,text='Close Window',command=_quit_sub)
            close_button.grid(row=8,column=2,columnspan=2)
            cwin.mainloop()
            
            
        
        custmap_button = Tk.Button(master=self.root,text = 'Custom Map', command = _customize)
        custmap_button.grid(row=5,column=1)
        quit_button = Tk.Button(master=self.root,text="Quit",command=_quit_win)
        quit_button.grid(row=5,column=8)



        Tk.mainloop()
        
        
def xyz(X,Y,Z):
     xnow = sum([complex(xi) for xi in X.split('+')])
     ynow = sum([complex(yi) for yi in Y.split('+')])
     znow = sum([complex(zi) for zi in Z.split('+')])
     xr,xi,yr,yi,zr,zi = float(np.real(xnow)),float(np.imag(xnow)),float(np.real(ynow)),float(np.imag(ynow)),float(np.real(znow)),float(np.imag(znow))
     return np.array([xr+1.0j*xi,yr+1.0j*yi,zr+1.0j*zi])/np.linalg.norm(np.array([xr+1.0j*xi,yr+1.0j*yi,zr+1.0j*zi]))
    