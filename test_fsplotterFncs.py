# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 12:46:46 2022

@author: Andrew

Idea is to allow testing of wien2k_lib.py functions
"""

import numpy as np
import chinook.build_lib as build_lib
import chinook.operator_library as operator_library
from chinook.ARPES_lib import experiment

import chinook.wien2k_lib as wien
import chinook.wien2k_plotting as wienPlot


import chinook.atomic_mass as am

import matplotlib.pyplot as plt



#case = 'data/FeTe_19nov21/FeTe_19nov21'
case = 'data/FeSe_mp20311/FeSe_mp20311'


avec,a,b,c = wien.get_PLV_from_struct(case)
k_object = wien.create_kobject(case)

# load in QTL, comment out if not nessecary as very memory intensive
QTL,bands,orbitals,E_F = wien.load_qtl(case)

'''
# FIRST TEST:
    # are the lengths of the kpath and bands data the same?
if len(k_object.kcut) == len(bands[1]):
    print('Test 1: kpath and band have same number of points.')
'''

'''
# SECOND TEST:
    # can we plot the bands on a figure, against k-path
f = plt.figure(0,(10,10))
ax = f.add_subplot()
for b in bands.keys():
    ax.plot(k_object.kcut,bands[b])

ax.set_ylim([E_F - 6, E_F+4])
ax.axhline(E_F,linestyle='--')
plt.show()
'''

'''

# THIRD TEST:
    # scatter plot with a qtl as a linewidth to see if that works...
    
w = lambda x: np.power(x,0.8)
    
f = plt.figure(0,(10,10))
ax = f.add_subplot()
for b in bands.keys():
    ax.scatter(k_object.kcut,bands[b],w(QTL[b][1][8]))#,alpha=QTL[b][1][8]) # 0 is 'tot' band

ax.set_ylim([E_F - 5, E_F+5])
ax.axhline(E_F,linestyle='--')
plt.show()
'''

'''
# Fourth Test:
    # using the wien2k_plotting plot_bands function
    
wienPlot.plot_bands(k_object, bands,E_F)
'''

'''
# FIFTH TEST:
    # using wien2k_plotting plot_bandCharacter function
orb = [1,7]
wienPlot.plot_bandCharacter(k_object, bands,QTL,orb)
wienPlot.plot_bandCharacter(k_object, bands, QTL, orb,w=lambda x: np.power(x,2))
wienPlot.plot_bandCharacter(k_object, bands, QTL, orb,w=lambda x: np.power(x,0.2))
'''

#SIXTH TEST:
    # using the wien2k_plotting plot_dominant_bandCharacter() function
proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
w = lambda x: 100*np.power(x,0.3)/np.max(x)
wienPlot.plot_dominant_bandCharacter(k_object, bands, QTL, proj,E_F=E_F)#,w=w)