#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 12:17:01 2017

@author: ryanday
"""
import numpy as np
import matplotlib.pyplot as plt
import build_lib
import ARPES_lib as ARPES
import scipy.ndimage as nd
import operator_library as ops
import Tk_plot

if __name__=="__main__":
    a,c =  3.7734,5.5258 
    avec = np.array([[a/np.sqrt(2),a/np.sqrt(2),0.0],[-a/np.sqrt(2),a/np.sqrt(2),0.0],[0.0,0.0,c]])

    filenm = 'FeSe_RD.txt'
    CUT,REN,OFF,TOL=a*3,1,0.15,0.001
    G,X,M,Z,mM = np.array([0,0,0]),np.array([0,0.5,0]),np.array([0.5,0.5,0]),np.array([0,0,0.5]),np.array([-0.5,-0.5,0])

	

    spin = {'soc':True,'lam':{0:0.04}}
    
    slab_dict = {'bool':False,
                'hkl':np.array([0,0,1]),
                'cells':20,
                'buff':8,
                'term':0,
                'avec':avec}

    Bd = {'atoms':[0,0],
			'Z':{0:26},
			'orbs':[["32xy","32XY","32xz","32yz","32ZR"],["32xy","32XY","32xz","32yz","32ZR"]],
			'pos':[np.array([-a/np.sqrt(8),0,0]),np.array([a/np.sqrt(8),0,0])],
            'slab':slab_dict}

    Kd = {'type':'F',
			'pts':[M,G,M],
			'grain':200,
			'labels':['M','$\Gamma$','M']}


    Hd = {'type':'txt',
			'filename':filenm,
			'cutoff':CUT,
			'renorm':REN,
			'offset':OFF,
			'tol':TOL,
			'so':spin['soc']}
 
    	#####
    Bd = build_lib.gen_basis(Bd,spin)
    Kobj = build_lib.gen_K(Kd,avec)
    TB = build_lib.gen_TB(Bd,Hd,Kobj)
    TB.solve_H()
    TB.plotting(-1.5,0.5)
#    O = ops.LdotS(TB,axis=None,vlims=(-0.5,0.5),Elims=(-0.5,0.5))
##    
##    
    ARPES_dict={'cube':{'X':[-0.3,0.3,60],'Y':[-0.3,0.3,60],'kz':0.0,'E':[-0.3,0.05,120]},
                'SE':[0.005,0.01],
                'directory':'C:\\Users\\rday\\Documents\\TB_ARPES\\2018\\TB_ARPES_2018\\FeSe',
                'hv': 21.2,
                'pol':np.array([0,np.sqrt(0.5),-np.sqrt(0.5)]),
                'mfp':7.0,
                'resolution':{'E':0.03,'k':0.05},
                'T':[True,10.0],
                'W':4.0,
                'angle':0.0,
                'spin':None,
                'slice':[False,-0.2]}

#
#    

    expmt = ARPES.experiment(TB,ARPES_dict)
    expmt.datacube(ARPES_dict)
    
    expmt.plot_gui(ARPES_dict)


    
    
	#####