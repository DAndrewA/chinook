# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 16:20:49 2022

@author: Andrew
"""

import copy

import numpy as np
import matplotlib.pyplot as plt
import distinctipy # used for distinct colours for bands
from matplotlib.lines import Line2D # used for custom legend

import chinook.wien2k_lib as wien




def plot_bands(k_obj,bands,E_F=0,ax=None,window=[-6,4]):
    '''
    Will plot all the bands for the bands object on a single axis
    '''
    if not ax: # if no axis supplied, create new set to plot on
        fig = plt.figure(0,(10,10))
        ax = fig.add_subplot()
    
    for band in bands.keys():
        ax.plot(k_obj.kcut,bands[band])
        
    ax.set_ylim([E_F + window[0], E_F + window[1]])
    ax.axhline(E_F,linestyle='--')
    
    ax = label_axis(ax, k_obj)
    #plt.show()
    return


def plot_bandCharacter(k_obj,bands,QTL,orb,E_F=None,ax=None,w= lambda x: x,window=[-6,4]):
    '''
    Creates a scatter plot for the bands with the linewidth being controlled by the scaling function w
    
    orb [1x2 array]: [atom,orbital] to be plotted (as indexed in QTL dictionary)
    '''
    if not ax:
        fig = plt.figure(0,(10,10))
        ax=fig.add_subplot()
        
    if not E_F:
        E_F = 0
    
    for b in bands.keys():
        ax.scatter(k_obj.kcut,bands[b],w(QTL[b][orb[0]][orb[1]]))
        
    
    ax.set_ylim([E_F + window[0], E_F+ window[1]])
    ax.axhline(E_F,linestyle='--')
    
    ax = label_axis(ax, k_obj)
    #plt.show()
    return
        
    
        



def plot_dominant_bandCharacter(k_obj,bands,QTL,proj,E_F=None,ax=None,w = lambda x: 1,window=[-6,4],colours = None,suppressLegend=False):
    '''
    Function to plot the dominant band character for a kpath
    
    k_obj: k-object as used in the chinook library
    bands: bands dictionary loaded in from load_qtl()
    QTL: the partial charges dictionary as loaded in from load_qtl()
    
    proj: array of lists, each element is an orbital to be considered.
    
    i.e. proj = [ [1,7,'$D_{xy}$'] , [1,6,'$D_{z^2}$'] ]
    would give the projection to be the 1st atom's 8th orbital, labelled as D_xy, etc
    '''
    if not ax: # if no axis supplied, create new figure to plot on
        fig = plt.figure(0,(10,10))
        ax = fig.add_subplot()
        
    if not E_F:
        E_F = 0
       
    # use the distinctipy package to generate unique colours for the orbitals
    if not colours:
        colours = distinctipy.get_colors(len(proj))
    if len(colours) < len(proj):
        colours = distinctipy.get_colors(len(proj))
       
    for b in bands.keys():
        # for each band, we need to calculate the dominant orbital character
        maxQTL = QTL[b][proj[0][0]][proj[0][1]]
        for p in proj:
            atom,orb = p[0],p[1]
            # for each projected orbital, finds the indices where its the dominant orbital and puts those values into maxQTL
            index_higher_QTL = QTL[b][atom][orb] > maxQTL
            maxQTL[index_higher_QTL] = QTL[b][atom][orb][index_higher_QTL]
            
        for i,p in enumerate(proj):
            atom,orb = p[0],p[1]
            index_dominant = QTL[b][atom][orb] >= maxQTL
            
            kpos = k_obj.kcut[index_dominant]
            E = bands[b][index_dominant]
            Q = QTL[b][atom][orb][index_dominant]
            ax.scatter(kpos,E,w(Q),color=colours[i])
            
        if b%10 == 0:print('band {} drawn'.format(b))
            
    
    ax.set_ylim([E_F + window[0], E_F+window[1]])
    ax.axhline(E_F,linestyle='--')
    
    if not suppressLegend:
        #create the custom legend for each of the orbitals:
        legend_elements = []
        for i,p in enumerate(proj):
            elem = Line2D([0],[0],color=colours[i],lw=5,label=p[2])
            legend_elements.append(elem)
        ax.legend(handles=legend_elements,loc='upper right')
    
    ax = label_axis(ax,k_obj)
    #plt.show()
    return ax
    
            
        
def plot_adjacent_dominant_bandCharacter(centre,direction,k_obj,bands,QTL,proj,E_F=None,ax=None,w = lambda x: 1,window=[-2,0.5],verbose=0,colours = None,suppressLegend=False):
    '''
    This function is to allow the easy plotting of a symmetric plot centred around a symmetry point, going in a specific direction
    '''    
    
    new_bands = {}
    new_QTL = {}
    ''' # SECTION NEEDS FIXING, CURRENTLY SELECTING ALL BANDS
    for b in bands.keys():
        maxE,minE = np.max(bands[b]),np.min(bands[b])
        if minE <= E_F + window[1] or maxE >= E_F + window[1]: # if a band comes within the allowed region of the Fermi energy
            new_bands[b] = bands[b]
    '''
    new_bands = copy.copy(bands)

    #print('accepted bands: {}'.format(new_bands.keys()))
    # get the indices for info between gamma and M
    labels = k_obj.labels
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
            ind_centre_k = list(k_obj.kcut).index( k_obj.kcut_brk[ind_adjacent[0]] )
            ind_direction_k = list(k_obj.kcut).index( k_obj.kcut_brk[ind_adjacent[1]] )
            
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
            new_k_object = copy.copy(k_obj)
            new_k_object.labels = [ labels[ind_adjacent[1]] , labels[ind_adjacent[0]], labels[ind_adjacent[1]] ]
            new_k_object.kcut = k_obj.kcut[index_band] - k_obj.kcut[index_band[0]] # centres k-values on 0, so we can have it go from -k_dir to +k_dir
            new_k_object.kcut = np.hstack(  ((new_k_object.kcut[-1:-len(new_k_object.kcut):-1] * -1), new_k_object.kcut)   )
            new_k_object.kcut_brk = [new_k_object.kcut[0]  ,  new_k_object.kcut[len(index_band)] ,new_k_object.kcut[-1]]

            if verbose: print('len(index_band2 = {})'.format(len(index_band2)))
            if verbose: print('len(new_k_object.kcut) = {}'.format(len(new_k_object.kcut)))

            # should be ready to plot!
                     
            #proj = [[1,6,'$D_{z^2}$'],[1,7,'$D_{xy}$'],[1,8,'$D_{x^2 y^2}$']]#,[1,9,'$D_{xz}+D_{yz}$']]
            #w = lambda x: 100*np.power(x,0.3)/np.max(x)
            ax = plot_dominant_bandCharacter(new_k_object, new_bands, new_QTL, proj,E_F=E_F,window=window,w=w,ax=ax, colours = colours,suppressLegend=suppressLegend)
            return ax
            
        else:
            print('No adjacent k-points {} and {} found in band'.format(centre,direction))
            return ax
    else:
        print('Specified points {}, {} not in k_object.labels'.format(centre,direction))
        return ax
    
        
   
    
   
    
   
def label_axis(ax,k_object):
    '''
    function to apply the axis labels from the k_object to the current axis
    '''
    labels = k_object.labels 
    # sets the x-tiks and labels appropriately
    ax.set_xticks(k_object.kcut_brk)
    ax.set_xticklabels(labels)
    for k in k_object.kcut_brk:
        ax.axvline(k)   
    return ax