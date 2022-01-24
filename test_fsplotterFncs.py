# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 12:46:46 2022

@author: Andrew

Idea is to allow testing of fsplotter_lib.py functions
"""

import numpy as np
import chinook.build_lib as build_lib
import chinook.operator_library as operator_library
from chinook.ARPES_lib import experiment

import chinook.fsplotter_lib as fs


import chinook.atomic_mass as am

import matplotlib.pyplot as plt



case = 'data/FeTe_19nov21'


avec,a,b,c = fs.get_PLV_from_struct(case)
k_object = fs.create_kobject(case)

# load in QTL, comment out if not nessecary as very memory intensive
QTL,bands,orbitals,E_F = fs.load_qtl(case)

'''
# FIRST TEST:
    # are the lengths of the kpath and bands data the same?
if len(k_object.kcut) == len(bands[1]):
    print('Test 1: kpath and band have same number of points.')
'''

'''
# SECOND TEST:
    # can we plot the bands on a figure, against k-path

ax = plt.axes()
for b in bands.keys():
    ax.plot(k_object.kcut,bands[b])

ax.set_ylim([E_F - 5, E_F+5])
ax.axhline(E_F,linestyle='--')
plt.show()

'''

'''
# THIRD TEST:
    # scatter plot with a qtl as a linewidth to see if that works...
    
w = lambda x: np.power(x,0.8)
    
ax = plt.axes()
for b in bands.keys():
    ax.scatter(k_object.kcut,bands[b],w(QTL[b][1][7]))#,alpha=QTL[b][1][8]) # 0 is 'tot' band

ax.set_ylim([E_F - 5, E_F+5])
ax.axhline(E_F,linestyle='--')
plt.show()
'''
# plot bands from bands and k_object data:
