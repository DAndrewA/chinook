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
print('a,b,c = {},{},{}'.format(a,b,c))
k_object = wien.create_kobject(case)

# load in QTL, comment out if not nessecary as very memory intensive
QTL,bands,orbitals,E_F = wien.load_qtl(case,verbose=0)


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

orbs_Fe = ['32ZR','32XY', '32xy','32xz', '32yz']
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
                    'tol':-1} # 1e-4
TB = build_lib.gen_TB(basis,hamiltonian_args,k_object)


'''
Having created the TB model and loaded in all the energy data, we want to see 
if we can get the bands and QTLs into the Eband and Evec data objects for the 
TB data object.
'''
TB.solve_H() # remember this is what we want to get around...

# this is a subset of bands that crosses the fermi level
bands_subset = [22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37]
band0 = bands_subset[0] # bands_subset needs to be ordered and consecutive

# put the bands into the TB object
shape_eband = (len(k_object.kpts),len(bands_subset))
Eband = np.zeros(shape_eband)
for b in bands_subset:
    print('Inserting band {} as band {}'.format(b,b-band0))
    for j,k in enumerate(k_object.kpts):
        Eband[j][b-band0] = bands[b][j]
TB.Eband = Eband

# put the QTLs into the Evec object
shape_evec = (len(k_object.kpts),len(bands_subset),len(TB.basis))
Evec = np.zeros(shape_evec)

# this will be an array simillar to proj from test_fsplotterFncs, in that it
# will show the way to get the desired orbital characters from the QTL dictionary's
# arrays.
wienOrbs = [[1,7], [1,8], [1,9], [1,10], [1,11], # the iron orbitals in the QTL dictionary
            [1,7], [1,8], [1,9], [1,10], [1,11], # duplicated twice because of two Fe in unit cell
            [2,3],[2,4],[2,5], # same for selenium, considering px py and pz
            [2,3],[2,4],[2,5]]

w = lambda x: x/2#np.power(x/2,2) # x/2
for b in bands_subset:
    for j,k in enumerate(k_object.kpts):
        for i,proj in enumerate(wienOrbs):
            Evec[j][b-band0][i] = w(QTL[b][proj[0]][proj[1]][j])
TB.Evec = Evec


TB.plotting(win_min=E_F-8,win_max=E_F+5) # the energy scale needed to see the first few bands...

dz2 = operator_library.fatbs(proj=[[0,5]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)
dx2y2 = operator_library.fatbs(proj=[[1,6]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)
dxy = operator_library.fatbs(proj=[[2,7]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)
dxz = operator_library.fatbs(proj=[[3,8]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)
dyz = operator_library.fatbs(proj=[[4,9]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)

px = operator_library.fatbs(proj=[[10,13]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)
py = operator_library.fatbs(proj=[[11,14]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)
pz = operator_library.fatbs(proj=[[12,15]],TB=TB,Elims=(E_F-8,E_F+5),degen=True)



'''
For the ARPES experiment simulation, we need to define a few variables before we can even begin:
    
    'cube': dictionary containing the region in k-space in which we want to view the bands.
            Distances are given in 1/Ang, and the final number is the number of points in that direction to calculate at.
            
(the rest are mostly self explanatory)

Our aim is to edit the actual ARPES library files to see how we can avoid having 
to calculate the band structure again in the region...
'''

arpes = {'cube':{'X':[-0.628,0.628,300],'Y':[-0.628,0.628,300],'E':[-0.05,0.05,50],'kz':0.0}, #domain of interest
    'hv':100,                          #photon energy, in eV
    'T':10,                            #temperature, in K
    'W':E_F,                           # work function, in this case I've used fermi energy
    'pol':np.array([1,0,-1]),           #polarization vector
    'SE':['constant',0.02],            #self-energy, assume for now to be a constant 20 meV for simplicity
    'resolution':{'E':0.02,'k':0.02}}  #resolution

arpes_experiment = experiment(TB,arpes) #initialize experiment object
#arpes_experiment.datacube() #execute calculation of matrix elements

#I,Ig,ax = arpes_experiment.spectral(slice_select=('w',0.0))