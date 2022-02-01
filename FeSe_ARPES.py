# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 16:19:49 2022

@author: Andrew


This is the code in which we'll aim to create our ARPES plots for the given case
"""

# firstly, we'll load in all the lattice parameters and objects required to create the band structure, etc.

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
#QTL,bands,orbitals,E_F = wien.load_qtl(case,verbose=0)


'''
In the basis object, we need the orbitals for each atom.

In basis['atoms'], we have an array of numbers indicating the species label for a particular atom
e.g. basis['atoms'] = [species of atom 0, species of atom 1,...]
                    = [0,0,1,1]
                    (implies the first two atoms are species 0, and the next two are for species 1)
                    

In basis['Z'], we have a dictionary describing the Z values for the atomic species
e.g. basis['Z'] = {species0:Z0 , species1:Z1, ...}
                = {0:26,1:34}
                (the first species is iron, and the next is selenium)
                
orbs: We want a full set of d orbitals for the iron atoms, and the p orbitals for the Se atoms
thus, let:
    orbs_Fe = ['32xy','32yz','32xz','32ZR','32XY']
    orbs_Se = ['41x','41y','41z']
    basis['orbs'] = [orbs_Fe,orbs_Fe,orbs_Se,orbs_Se]
    
In basis['pos'], we have an array of numpy arrays that contain the positions of the atoms in Cartesian coordinates in units of angstrom
Positions should be extracted from the .struct file already.


SPIN ORBIT COUPLING:
    For now, we'll focus on the non-SO case. However, if we want to include SO coupling, will need to determine a coupling factor lambda for each atomic species.
 
spin = {'bool':False}
#spin = {'bool':True,  #include spin-degree of freedom: double the orbital basis
#'soc':True,    #include atomic spin-orbit coupling in the calculation of the Hamiltonian
#'lam':{0:0.5}} #spin-orbit coupling strength in eV, require a value for each unique species in our basis

'''
# creating the basis object for the Chinook data

orbitalNames_wien = ['0','PX','PY','PZ','DZ2','DX2Y2','DXY','DXZ','DYZ']
orbitalNames_chinook = ['0','1x','1y','1z','2ZR','2XY','2xy','2xz','2yz']

orbs_Fe = ['32xy','32yz','32xz','32ZR','32XY']
orbs_Se = ['41x','41y','41z']
orbs = [orbs_Fe,orbs_Fe,orbs_Se,orbs_Se]

basis,spin = wien.create_basisObject(case, avec,orbs=orbs)


'''
For the TB model, we don't care how the potential/hopping terms are set up, as long
as we can avoid using them later. As such, as a first effort, I'll work using a
nearest neighbour model (simplest) and populate every hopping term with the number 0.

In future, we can test if ew are avoiding the use of these terms by populating
the model with random energies and seeing if we change the ARPES plot.
'''
# creating the TB model for our system

V_SK = {}
for i,a in enumerate(basis['orbs']):
    for o in a:
        # for each orbital for each atom, define the on-site energy as 0
        orbitalString = '{}{}'.format(i,o[0:2]) # take first two characters as want it for 4p in general, or 3d, etc
        V_SK[orbitalString] = 0

hamiltonian_args = {'type':'SK',
                    'V':V_SK,
                    'avec':avec,
                    'cutoff':0, # 2.5
                    'renorm':0, # 1.0
                    'offset':0.0,
                    'tol':1e-3} # 1e-4
TB = build_lib.gen_TB(basis,hamiltonian_args)
    
