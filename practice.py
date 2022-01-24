# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:14:01 2022

@author: Andrew
"""

import numpy as np
import chinook.build_lib as build_lib
import chinook.operator_library as operator_library
from chinook.ARPES_lib import experiment

import chinook.wien2k_lib as wien


import chinook.atomic_mass as am



case = 'data/FeTe_19nov21'



avec,a,b,c = wien.get_PLV_from_struct(case)



#--------------------- DEFNING kpath ------------------------------------------


'''
# Define the k path. In this case, we'll go gamma -> Z -> corner -> gamma
kpoints = np.array([[0,0,0],[0,0,0.5],[0.5,0.5,0.5],[0,0,0]])
labels = np.array(['$\\Gamma$','$Z$','corner','$\\Gamma$'])

kdict = {'type':'F', #fractional coordinates for kpoints
'avec':avec,
'pts':kpoints, #pts is the key points along our path
'grain':200, # number of points along each part of path?
'labels':labels}

k_object = build_lib.gen_K(kdict)
'''
k_object = wien.create_kobject(case)


#---------------------------- DEFINING STRUCTURE -----------------------------
'''
spin = {'bool':True,  #include spin-degree of freedom: double the orbital basis
'soc':True,    #include atomic spin-orbit coupling in the calculation of the Hamiltonian
'lam':{0:0.5}} #spin-orbit coupling strength in eV, require a value for each unique species in our basis


Sb1 = np.array([0.0,0.0,0.0])
Sb2 = np.array([np.sqrt(0.5)*a,0,0])


basis = {'atoms':[0,0], #two equivalent atoms in the basis, both labelled as species #0
 'Z':{0:51},     #We only have one atomic species, which is antimony #51 in the periodic table.
 'orbs':[['51x','51y','51z'],['51x','51y','51z']], #each atom includes a full 5p basis in this model, written in n-l-xx format
 'pos':[Sb1,Sb2], #positions of the atoms, in units of Angstrom
 'spin':spin} #spin arguments.

basis_object = build_lib.gen_basis(basis)
'''
basis_object,spin = wien.create_basisObject(case, avec)


#------------------------------------------------------------------------------


# DEFININF THE TB HAMILTONIAN, THIS IS THE BIT WE WANT TO SKIP
Ep = 0.7
Vpps = 0.25
Vppp = -1.0
V1 = {'051':Ep,'005511S':Vpps,'005511P':Vppp}
V2 = {'005511S':Vpps/a,'005511P':Vppp/a}
VSK = [V1,V2]

cutoff = [0.8*c,1.1*c]

hamiltonian = {'type':'SK',     #Slater-Koster type Hamiltonian
      'V':VSK,          #dictionary (or list of dictionaries) of onsite and hopping potentials
       'avec':avec,     #lattice geometry
      'cutoff':cutoff,  #(cutoff or list of cutoffs) maximum length-scale for each hoppings specified by VSK
      'renorm':1.0,     #renormalize bandwidth of Hamiltonian
       'offset':0.0,    #offset the Fermi level
      'tol':1e-4,       #minimum amplitude for matrix element to be included in model.
      'spin':spin}      #spin arguments, as defined above


# Ultimately, we need the TB object, but we want to automate its cration and remove as much calculation as possible.
TB = build_lib.gen_TB(basis_object,hamiltonian,k_object)
# could use case.struct to generate basis object...


#----------------- THIS IS THE BIT WHERE WE'LL SNEAKY OUR BAND STRUCTURE IN -----
TB.Kobj = k_object
TB.solve_H()


# need to figure out how to get correct bands format into TB data structure

#k_object,TB = fs.change_kobject_TB(case, k_object, TB)

TB.plotting(win_min=-14.0,win_max=8.0)



#--------------- THIS IS THE BIT WHERE THE MAGIC HAPPENS!!! ----------
arpes = {'cube':{'X':[-0.628,0.628,300],'Y':[-0.628,0.628,300],'E':[-0.05,0.05,50],'kz':0.0}, #domain of interest
    'hv':100,                          #photon energy, in eV
     'T':10,                           #temperature, in K
    'pol':np.array([1,0,-1]),           #polarization vector
    'SE':['constant',0.02],            #self-energy, assume for now to be a constant 20 meV for simplicity
    'resolution':{'E':0.02,'k':0.02}}  #resolution

arpes_experiment = experiment(TB,arpes) #initialize experiment object
arpes_experiment.datacube() #execute calculation of matrix elements

I,Ig,ax = arpes_experiment.spectral(slice_select=('w',0.0))