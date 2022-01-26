# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 12:46:46 2022

@author: Andrew

Idea is to allow testing of wien2k_lib.py functions
"""
import distinctipy as dispy

import numpy as np
import chinook.build_lib as build_lib
import chinook.operator_library as operator_library
from chinook.ARPES_lib import experiment

import chinook.wien2k_lib as wien
import chinook.wien2k_plotting as wienPlot

from matplotlib.lines import Line2D # used for custom legend

import chinook.atomic_mass as am

import matplotlib.pyplot as plt



#case = 'data/FeTe_19nov21/FeTe_19nov21'
#case = 'data/FeSe_mp20311/FeSe_mp20311'
case = 'data/FeSe_25jan22/FeSe_25jan22'


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

'''
#SIXTH TEST:
    # using the wien2k_plotting plot_dominant_bandCharacter() function
#bands[59] = 2*(bands[59]-E_F) + E_F   

proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
w = lambda x: 100*np.power(x,0.3)/np.max(x)
wienPlot.plot_dominant_bandCharacter(k_object, bands, QTL, proj,E_F=E_F)#,w=w)
'''

'''
# SEVENTH TEST:
    # some basic band shifting, will raise lower every odd band by 0.1ev, raise every even band by 0.1ev
for b in bands.keys():
    if b%2 == 0:
        bands[b] = bands[b] + 0.02
    else:
        bands[b] = bands[b] - 0.02

proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
w = lambda x: 100*np.power(x,0.3)/np.max(x)
wienPlot.plot_dominant_bandCharacter(k_object, bands, QTL, proj,E_F=E_F)#,w=w)
'''
'''
# EIGTH TEST:
    # Combine each doubling of bands, as well as their QTL's into a new bands and QTL data structure
new_bands = {}
new_QTL = {}
for b in bands.keys():
    if b%2 == 1: # if the first band in a pair
        temp_band = bands[b]
        temp_QTL = QTL[b]
    else: # on the second band in the pair
        new_bands[int(b/2)] = temp_band
        for atom in temp_QTL.keys():
            temp_QTL[atom] = (temp_QTL[atom] + QTL[b][atom])/2 # average of two QTLs
            new_QTL[int(b/2)] = temp_QTL
        print('bands {} and {} combined'.format(b-1,b))
            
new_bands[29] = 2*(new_bands[29] - new_bands[29][0]) + new_bands[29][0]
        
proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
w = lambda x: 100*np.power(x,0.3)/np.max(x)
wienPlot.plot_dominant_bandCharacter(k_object, new_bands, new_QTL, proj,E_F=E_F)#,w=w)
'''  

'''
# NINTH TEST:
    # going to plot M-GAMMA-M within +1,-2ev of the fermi level

new_bands = {}
new_QTL = {}
window = [-2,0.5]
for b in bands.keys():
    maxE,minE = np.max(bands[b]),np.min(bands[b])
    if minE <= E_F + window[1] or maxE >= E_F + window[1]: # if a band comes within the allowed region of the Fermi energy
        new_bands[b] = bands[b]
        new_QTL[b] = QTL[b]
        

#print('accepted bands: {}'.format(new_bands.keys()))
# get the indices for info between gamma and M
centre = '$Z$'
direction = '$A$'
labels = k_object.labels
if centre in labels and direction in labels: # need to check the desired direction is include in the band
    # get the indices of the centre and direction. If two are adjacent in the array, there's a direct calculation between them
    ind_centre = [i for i in range(len(labels)) if labels[i] == centre]
    ind_direction = [i for i in range(len(labels)) if labels[i] == direction]
    
    ind_adjacent = []
    # need to check for adjacency in these indexes
    for i in ind_centre:
        for j in ind_direction:
            if i-j == 1 or j-i == 1: # if adjacent
                ind_adjacent = [i,j]
                break
        if ind_adjacent: # if indices have been found, break from searching
            break
    
    # only if adjacent indices have been found:
    if ind_adjacent:
        print('adjacent indices found: [{},{}]'.format(ind_adjacent[0],ind_adjacent[1]))
        # select the indices for the bands in which the information lies
        ind_centre_k = list(k_object.kcut).index( k_object.kcut_brk[ind_adjacent[0]] )
        ind_direction_k = list(k_object.kcut).index( k_object.kcut_brk[ind_adjacent[1]] )
        
        index_band = []
        if ind_adjacent[0] < ind_adjacent[1]:
            index_band = list(range(ind_centre_k,ind_direction_k + 1))
        else:
            index_band = list(range(ind_centre_k,ind_direction_k - 1, -1))
            
        #print('index_band = {}'.format(index_band))
        # the current index list goes from the index for the centre to that of the 'direction'
        # we want to reverse this and prepend in to the set of indices
        index_band2 = index_band[-1:-len(index_band):-1]
        index_band2 = np.hstack( (index_band2,index_band) )
        #print('index_band2 step1 = {}'.format(index_band[-1:-len(index_band):-1]))
        
        #print('index_band2 = {}'.format(index_band2))
        #index_band2 = np.array(index_band2) # turn into numpy array for convenience
        index_band = np.array(index_band) # distinction made between doubled up band and original fo use in extracting k values later
        
        #print('list(index_band2) = {}'.format(list(index_band2)))
        # goes through each band in new_bands and selects elements within our range
        '''
'''
        for b in new_bands.keys():
            new_bands[b] = new_bands[b][index_band2]
            print('new band:')
            print('len(new_bands[{}]) = {}'.format(b,len(new_bands[b])))
            print('len(bands[{}]) = {}'.format(b,len(bands[b])))
            for atom in new_QTL[b].keys(): # also selects the QTL elements in desired range
                # need to create new matrix for qtls of size (N orbitals)x(len(index_band2))
                shape_temp = ( np.shape(new_QTL[b][atom])[0] , len(index_band2) )
                #print('new atom, shape_temp = {}'.format(shape_temp))
                print('new atom, shape(new_QTL[{}][{}]) = {}'.format(b,atom,np.shape(new_QTL[b][atom])))
                print('new atom, shape(QTL[{}][{}]) = {}'.format(b,atom,np.shape(QTL[b][atom])))
                new_QTL_temp = np.zeros(shape_temp)
                for orb in range( np.shape(new_QTL[b][atom])[0] ):
                    print('len(new_QTL_temp[{}]) = {}'.format(orb,len(new_QTL_temp[orb])))
                    print('len(new_QTL[{}][{}][{}]) = {}'.format(b,atom,orb,len(new_QTL[b][atom][orb])))
                    print('len(new_QTL[{}][{}][{}][index_band2]) = {}'.format(b,atom,orb,len(new_QTL[b][atom][orb][index_band2])))
                    new_QTL_temp[orb] = new_QTL[b][atom][orb][index_band2]
                new_QTL[b][atom] = new_QTL_temp
'''
'''
        #-------  DOING BANDS AND QTL SEPERATELY, FOR NOW
        for b in new_bands.keys():
            #l1 = len(new_bands[b])
            new_bands[b] = new_bands[b][index_band2]
            #l2 = len(new_bands[b])
            #print('Band reshaping, band {}, {} -> {}'.format(b,l1,l2))
            
        new_QTL = {}
        for b in new_bands.keys():
            band_QTL = {}
            for atom in QTL[b].keys(): # for each atom in the original band
                shape_temp = ( np.shape(QTL[b][atom])[0] , len(index_band2) )
                QTL_temp = np.zeros(shape_temp)
                for orb in range( np.shape(QTL[b][atom])[0] ):
                    QTL_temp[orb] = QTL[b][atom][orb][index_band2]
                band_QTL[atom] = QTL_temp
            new_QTL[b] = band_QTL
        
        # need to alter kobject to have appropriate k-values, labels, etc
        new_k_object = k_object
        new_k_object.labels = [ labels[ind_adjacent[1]] , labels[ind_adjacent[0]], labels[ind_adjacent[1]] ]
        new_k_object.kcut = k_object.kcut[index_band] - k_object.kcut[index_band[0]] # centres k-values on 0, so we can have it go from -k_dir to +k_dir
        new_k_object.kcut = np.hstack(  ((new_k_object.kcut[-1:-len(new_k_object.kcut):-1] * -1), new_k_object.kcut)   )
        new_k_object.kcut_brk = [new_k_object.kcut[0]  ,  new_k_object.kcut[len(index_band)] ,new_k_object.kcut[-1]]

        print('len(index_band2 = {})'.format(len(index_band2)))
        print('len(new_k_object.kcut) = {}'.format(len(new_k_object.kcut)))

        # should be ready to plot!
                 
        proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
        #w = lambda x: 100*np.power(x,0.3)/np.max(x)
        wienPlot.plot_dominant_bandCharacter(new_k_object, new_bands, new_QTL, proj,E_F=E_F,window=window)#,w=w)
        
        
    else:
        print('No adjacent k-points {} and {} found in band'.format(centre,direction))
'''

'''
#TENTH TEST;
    # should be ninth test, but as a function
proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
wienPlot.plot_adjacent_dominant_bandCharacter('$\Gamma$', '$M$', k_object, bands, QTL, proj,E_F=E_F)
'''

'''
# ELEVENTH TEST:
    # going to do a 4x4 subplot grid, each subplot containing the plot around a centre in a given direction

proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
colours = dispy.get_colors(len(proj))

# extracts unique labels
t = k_object.labels
print(t)
labels = []
for p in t:
    if p not in labels:
        labels.append(p)
plt.ion()
f = plt.figure(0,(13,13),600)
for x,centre in enumerate(labels):
    for y,direction in enumerate(labels):
        print('(x,y) = ({},{})'.format(x,y))
        print('centre = {} ; direction = {}'.format(centre,direction))
        print('k_object.labels = {}'.format(k_object.labels))
        ax = f.add_subplot(len(labels),len(labels),( len(labels)*x +y +1 ))
        ax = wienPlot.plot_adjacent_dominant_bandCharacter(centre, direction, k_object, bands, QTL, proj,E_F=E_F,ax=ax,colours=colours,suppressLegend=True)
        
        if x == 0 and y == 0: # create legend in upper right corner
            #create the custom legend for each of the orbitals:
            legend_elements = []
            for i,p in enumerate(proj):
                elem = Line2D([0],[0],color=colours[i],lw=5,label=p[2])
                legend_elements.append(elem)
            ax.legend(handles=legend_elements,loc='center')
plt.show()
'''