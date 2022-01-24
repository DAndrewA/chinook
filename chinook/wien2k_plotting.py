# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 16:20:49 2022

@author: Andrew
"""

import numpy as np
import matplotlib.pyplot as plt
import distinctipy # used for distinct colours for bands
from matplotlib.lines import Line2D # used for custom legend

import chinook.wien2k_lib as wien




def plot_bands(k_obj,bands,E_F=None,ax=None):
    '''
    Will plot all the bands for the bands object on a single axis
    '''
    if not ax: # if no axis supplied, create new set to plot on
        fig = plt.figure(0,(10,10))
        ax = fig.add_subplot()
    
    for band in bands.keys():
        ax.plot(k_obj.kcut,bands[band])
        
    if E_F: # if a fermi energy is supplied, centre the plot on it
        ax.set_ylim([E_F - 6, E_F+4])
        ax.axhline(E_F,linestyle='--')
    
    ax = label_axis(ax, k_obj)
    plt.show()
    return


def plot_bandCharacter(k_obj,bands,QTL,orb,E_F=None,ax=None,w= lambda x: x):
    '''
    Creates a scatter plot for the bands with the linewidth being controlled by the scaling function w
    
    orb [1x2 array]: [atom,orbital] to be plotted (as indexed in QTL dictionary)
    '''
    if not ax:
        fig = plt.figure(0,(10,10))
        ax=fig.add_subplot()
    
    for b in bands.keys():
        ax.scatter(k_obj.kcut,bands[b],w(QTL[b][orb[0]][orb[1]]))
        
    if E_F: # if a fermi energy is supplied, centre the plot on it
        ax.set_ylim([E_F - 6, E_F+4])
        ax.axhline(E_F,linestyle='--')
    
    ax = label_axis(ax, k_obj)
    plt.show()
    return
        
    
        



def plot_dominant_bandCharacter(k_obj,bands,QTL,proj,E_F=None,ax=None,w = lambda x: x):
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
       
    # use the distinctipy package to generate unique colours for the orbitals
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
            ax.scatter(kpos,E,w(Q),c=colours[i])
            
    #create the custom legend for each of the orbitals:
    legend_elements = []
    for i,p in enumerate(proj):
        elem = Line2D([0],[0],color=colours[i],lw=5,label=p[2])
        legend_elements.append(elem)
    ax.legend(handles=legend_elements)
    
    ax = label_axis(ax,k_obj)
    plt.show()
    return
    
            
        
   
    
   
    
   
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