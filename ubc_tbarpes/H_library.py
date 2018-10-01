#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 21:38:24 2017

@author: ryanday

Slater Koster Library for generating a Hamiltonian of Slater-Koster
MIT License

Copyright (c) 2018 Ryan Patrick Day

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import numpy as np
import ubc_tbarpes.orbital as orb
import ubc_tbarpes.SK as SK

hb = 6.626*10**-34/(2*np.pi)
c  = 3.0*10**8
q = 1.602*10**-19
A = 10.0**-10
me = 9.11*10**-31
mn = 1.67*10**-27
kb = 1.38*10**-23


#def slabify():
    

def txt_build(filename,cutoff,renorm,offset,tol):
    
    Hlist = []
    
    with open(filename,'r') as origin:
        
        for line in origin:
            
            spl = line.split(',')
            R = np.array([float(spl[2]),float(spl[3]),float(spl[4])])
            Hval = complex(spl[5])
            
            if len(spl)>7:
                Hval+=1.0j*float(spl[6])
            if abs(Hval)>tol and  np.linalg.norm(R)<cutoff:
                Hval*=renorm
                if np.linalg.norm(R)==0.0:
                    Hval-=offset
                    
                tmp = [int(spl[0]),int(spl[1]),R[0],R[1],R[2],Hval]

                Hlist.append(tmp)
            
    origin.close()
            
    return Hlist


                

def sk_build(avec,basis,V,cutoff,tol,renorm,offset,spin):
    '''
    Would like to find a better way of doing this, or at least getting around the whole cluster thing...
    '''
    try:
        reg_cut = [0.0]+list(cutoff)
        reg_cut = np.array(reg_cut)
    except TypeError:
        try:
            reg_cut = np.array([cutoff])
        except TypeError:
            print('Invalid cutoff-format')
            return None
    pt_max = np.ceil(np.array([(reg_cut).max()/np.linalg.norm(avec[i]) for i in range(len(avec))]).max())
    pts = region(int(pt_max)+1)
    H_raw = []
    o1o2norm = {}
    if spin:
        brange = int(len(basis)/2)
    else:
        brange = len(basis)
    for o1 in basis[:brange]:
        for o2 in basis[:brange]:
            if o1.index<=o2.index:
                for p in pts:
                    Rij = o2.pos-o1.pos+np.dot(p,avec)
            
                    Rijn = np.linalg.norm(Rij)
 

                    orb_label = "{:d}-{:d}-{:0.3f}-{:0.3f}-{:0.3f}".format(o1.index,o2.index,Rij[0],Rij[1],Rij[2])
                    if Rijn>max(reg_cut):

                        o1o2norm[orb_label] = True
                        mat_el=0.0
                
                    try:
                        o1o2norm[orb_label] #check to see if a pairing of these orbitals, along this direction, has already been calculated--if so, skip to next o1,o2                     
                        mat_el=0.0
                        continue
                    except KeyError: #if this pair has not yet been included, proceed

                        if isinstance(V,list): #if we have given a list of SK dictionaries (relevant to different distance ranges)                        
                            for i in range(len(reg_cut)-1):
                                if reg_cut[i]<=Rijn<reg_cut[i+1]: #if in this range of the dictionaries, use the lower bound
                                    tmp_V = V[i]
                                    mat_el = SK.SK_coeff(o1,o2,Rij,tmp_V,renorm,offset,tol)  #then matrix element is computed using the SK function
                        elif isinstance(V,dict): #if the SK matrix elements brought in NOT as a list of dictionaries...                               
                            mat_el = SK.SK_coeff(o1,o2,Rij,V,renorm,offset,tol)

                            
                    if abs(mat_el)>tol: 
                        H_raw.append([o1.index,o2.index,Rij[0],Rij[1],Rij[2],mat_el])                            
                    o1o2norm[orb_label] = True #now that the pair has been calculated, disqualify from subsequent calculations
    return H_raw


def sk_build_2(avec,basis,V,cutoff,tol,renorm,offset,spin):
    '''
    Testing a function to build SK model from using D-matrices, rather than the list of SK terms from table
    '''
    normal_order = {0:{'':0},1:{'x':0,'y':1,'z':2},2:{'xz':0,'yz':1,'xy':2,'ZR':3,'XY':4},3:{'z3':0,'xz2':1,'yz2':2,'xzy':3,'zXY':4,'xXY':5,'yXY':6}}
    reg_cut,pts = cluster_build(cutoff,avec)
    H_raw = []
    o1o2norm = {}
    if spin:
        brange = int(len(basis)/2)
    else:
        brange = len(basis)
    
    pair_channels = orbital_pair_channels(basis[:brange])
    SK_matrices = SK.SK_full(basis[:brange])
    for pc in pair_channels:
        
        Vmat = SK_matrices[pc[:4]]
            
        for p in pts:
            Rij = np.array(pc[-3:])+np.dot(p,avec)
            Rijn = np.linalg.norm(Rij)
            if Rijn<max(reg_cut):
                A,B = np.arctan2(Rij[1],Rij[0]),np.arccos(Rij[2]/Rijn)
                V = [0.1]*(2*(min(pc[1],pc[3]))+1)
                SKmat = Vmat(A,B,V)
                H_raw = H_raw + mat_els(Rij,SKmat)
                    
                
                
def mat_els(Rij,SKmat):
    return None
                    
    
            




def orbital_pair_channels(basis):
    '''
    Define the set of viable atom, orbital channels which can be connected via hopping. These determine the hopping matrices which need be defined.
    '''
    orb_chann = []
    for o in basis:
        if (o.atom,o.l) not in orb_chann:
            orb_chann.append((o.atom,o.l,*o.pos))
    pair_chann = [(orb_chann[yi][0],orb_chann[yi][1],orb_chann[yj][0],orb_chann[yj][1],*(orb_chann[yj][2:5]-orb_chann[yi][2:5])) for yi in range(len(orb_chann)) for yj in range(yi,len(orb_chann))]
    
    return pair_chann


def cluster_build(cutoff,avec):
    '''
    Generate a safe cluster of neighbouring lattice points to use in defining the hopping paths
    Return an array of lattice points which go safely to the edge of the cutoff range.
    
    '''
    try:
        reg_cut = [0.0]+list(cutoff)
        reg_cut = np.array(reg_cut)
    except TypeError:
        try:
            reg_cut = np.array([cutoff])
        except TypeError:
            print('Invalid cutoff-format')
            return None
    pt_max = np.ceil(np.array([(reg_cut).max()/np.linalg.norm(avec[i]) for i in range(len(avec))]).max())
    pts = region(int(pt_max)+1)
    return reg_cut,pts


def spin_double(H,lb):
    lenb = int(lb/2)
    h2 = []
    for i in range(len(H)):
        h2.append([H[i][0]+lenb,H[i][1]+lenb,H[i][2],H[i][3],H[i][4],H[i][5]])
    return h2


def SO(basis):
    '''
    Generate L.S  matrix-elements for a given basis. 
    This is generic to all l, except the normal_order, specifically defined for l=1,2...should generalize
    Otherwise, this structure holds for all l!
    
    In the factors dictionary, the weight of the different LiSi terms is defined. The keys are tuples of (L+/-/z,S+/-/z) in a bit
    of a cryptic way. For L, range (0,1,2) ->(-1,0,1) and for S range (-1,0,1) = S1-S2 with S1/2 = +/- 1 here
    
    L+,L-,Lz matrices are defined for each l shell in the basis, transformed into the basis of cubic harmonics.
    The nonzero terms will then just be used along with the spin and weighted by the factor value, and slotted into 
    a len(basis)xlen(basis) matrix HSO
    args:
        basis -- tight-binding model basis orbital list, as defined in TB_lib.py
    return:
        HSO list of matrix elements in standard format [i,j,0,0,0,HSO[i,j]]
    '''
    Md = Yproj(basis)
    normal_order = {0:{'':0},1:{'x':0,'y':1,'z':2},2:{'xz':0,'yz':1,'xy':2,'ZR':3,'XY':4},3:{'z3':0,'xz2':1,'yz2':2,'xzy':3,'zXY':4,'xXY':5,'yXY':6}}
    factors = {(2,-1):0.5,(0,1):0.5,(1,0):1.0}
    L,al = {},[]
    HSO = []
    for o in basis[:int(len(basis)/2)]:
        if (o.atom,o.l) not in al:
            al.append((o.atom,o.l))
            Mdn = Md[(o.atom,o.l,-1)]
            Mup = Md[(o.atom,o.l,1)]
            Mdnp = np.linalg.inv(Mdn)
            Mupp = np.linalg.inv(Mup)
            L[(o.atom,o.l)] = [np.dot(Mupp,np.dot(Lm(o.l),Mdn)),np.dot(Mdnp,np.dot(Lz(o.l),Mdn)),np.dot(Mdnp,np.dot(Lp(o.l),Mup))]

#        if (o.atom,o.l) not in al:
#            al.append((o.atom,o.l))
#            M = Md[(o.atom,o.l)]
#            Mp = np.linalg.inv(M)
#            L[(o.atom,o.l)] = [np.dot(Mp,np.dot(Lm(o.l),M)),np.dot(Mp,np.dot(Lz(o.l),M)),np.dot(Mp,np.dot(Lp(o.l),M))]

    for o1 in basis:
        for o2 in basis:
            if o1.index<=o2.index:
                LS_val = 0.0
                if np.linalg.norm(o1.pos-o2.pos)<0.0001 and o1.l==o2.l and o1.n==o2.n:
                    inds = (normal_order[o1.l][o1.label[2:]],normal_order[o2.l][o2.label[2:]])
                    
                    ds = (o1.spin-o2.spin)/2.
                    if ds==0:
                        s=0.5*np.sign(o1.spin)
                    else:
                        s=1.0
                    for f in factors:
                        if f[1]==ds:
                            LS_val+=o1.lam*factors[f]*L[(o1.atom,o1.l)][f[0]][inds]*s
                    HSO.append([o1.index,o2.index,0.,0.,0.,LS_val])

    return HSO


def Lp(l):
    '''
    L+ operator in the l,m_l basis, organized with (0,0) = |l,l>, (2*l,2*l) = |l,-l>
    The nonzero elements are on the upper diagonal
    arg: l int orbital angular momentum
    return M: numpy array (2*l+1,2*l+1) of real float
    '''
    M = np.zeros((2*l+1,2*l+1))
    r = np.arange(0,2*l,1)
    M[r,r+1]=1.0
    vals = [0]+[np.sqrt(l*(l+1)-(l-m)*(l-m+1)) for m in range(1,2*l+1)]
    M = M*vals
    return M

def Lm(l):
    '''
    L- operator in the l,m_l basis, organized with (0,0) = |l,l>, (2*l,2*l) = |l,-l>
    The nonzero elements are on the upper diagonal
    arg: l int orbital angular momentum
    return M: numpy array (2*l+1,2*l+1) of real float
    '''
    M = np.zeros((2*l+1,2*l+1))
    r = np.arange(1,2*l+1,1)
    M[r,r-1]=1.0
    vals = [np.sqrt(l*(l+1)-(l-m)*(l-m-1)) for m in range(0,2*l)]+[0]
    M = M*vals
    return M

def Lz(l):
    '''
    Lz operator in the l,m_l basis
    arg: l int orbital angular momentum
    retun numpy array (2*l+1,2*l+1)
    '''
    return np.identity(2*l+1)*np.array([l-m for m in range(2*l+1)])




def Yproj(basis):
    '''
    Define the unitary transformation rotating the basis of different inequivalent atoms in the
    basis to the basis of spherical harmonics for sake of defining L.S operator in basis of user
    args: basis--list of orbital objects
    
    returns: dictionary of matrices for the different atoms and l-shells--keys are tuples of (atom,l)
    
    Note this works only on p and d type orbitals, s is irrelevant, not currently supporting f orbitals
    
    29/09/2018 added reference to the spin character 'sp' to handle rotated systems effectively
    
    '''
    normal_order = {0:{'':0},1:{'x':0,'y':1,'z':2},2:{'xz':0,'yz':1,'xy':2,'ZR':3,'XY':4},3:{'z3':0,'xz2':1,'yz2':2,'xzy':3,'zXY':4,'xXY':5,'yXY':6}}
    a = basis[0].atom
    l = basis[0].l
    sp = basis[0].spin
    M = {}
    M_tmp = np.zeros((2*l+1,2*l+1),dtype=complex)
    for b in basis:
        label = b.label[2:]

        if b.atom==a and b.l==l and b.spin==sp:
            for p in b.proj:
                M_tmp[l-int(p[-1]),normal_order[l][label]] = p[0]+1.0j*p[1]
                
        else:
                #If we are using a reduced basis, fill in orthonormalized projections for other states in the shell
                #which have been ignored in our basis choice--these will still be relevant to the definition of the LS operator
            M_tmp = fillin(M_tmp,l)            
            M[(a,l,sp)] = M_tmp
                ##Initialize the next M matrix               
            a = b.atom
            l = b.l
            sp = b.spin
            M_tmp = np.zeros((2*l+1,2*l+1),dtype=complex)
            for p in b.proj:
                M_tmp[l-int(p[-1]),normal_order[l][label]] = p[0]+1.0j*p[1]
    
    M_tmp = fillin(M_tmp,l)
    M[(a,l,sp)] = M_tmp
    
    return M

def fillin(M,l):
    normal_order_rev = {0:{0:''},1:{0:'x',1:'y',2:'z'},2:{0:'xz',1:'yz',2:'xy',3:'ZR',4:'XY'},3:{0:'z3',1:'xz2',2:'yz2',3:'xzy',4:'zXY',5:'xXY',6:'yXY'}}

    for m in range(2*l+1):
        if np.linalg.norm(M[:,m])==0: #if column is empty (i.e. user-defined projection does not exist)
            proj = np.zeros(2*l+1,dtype=complex) 
            for pi in orb.projdict[str(l)+normal_order_rev[l][m]]: 
                proj[l-int(pi[-1])] = pi[0]+1.0j*pi[1] #fill the column with generic projection for this orbital (this will be a dummy)
            for mp in range(2*l+1): #Orthogonalize against the user-defined projections
                if np.linalg.norm(M[:,mp])!=0:
                    proj = GrahamSchmidt(proj,M[:,mp])
            M[:,m] = proj            
    return M



def AFM_order(basis,dS,p_up,p_dn):
  '''
  Add antiferromagnetism to the tight-binding model, by adding a different on-site energy to 
  orbitals of different spin character, on the designated sites. 
  args:
      basis -- list of orbital objects
      dS -- size of spin-splitting (eV) float
      p_up,p_dn -- numpy array of float indicating the orbital positions for the AFM order
  return: list of matrix elements, as conventionally arranged
  '''
  h_AF = []
  for bi in basis:
      if np.linalg.norm(bi.pos-p_up)==0:
          if bi.spin<0:
              h_AF.append([bi.index,bi.index,0,0,0,dS])
          else:
              h_AF.append([bi.index,bi.index,0,0,0,-dS])
      elif np.linalg.norm(bi.pos-p_dn)==0:
          if bi.spin<0:
              h_AF.append([bi.index,bi.index,0,0,0,-dS])
          else:
              h_AF.append([bi.index,bi.index,0,0,0,dS])
  return h_AF
    
    
def FM_order(basis,dS):
    '''
     Add ferromagnetism to the system. Take dS to assume that the splitting puts spin-up lower in energy by dS,
     and viceversa for spin-down. This directly modifies the TB_model's mat_els attribute
     args:
         basis -- list of orbital objects in basis
         dS -- float: energy of the spin splitting (eV)
    retun: list of matrix elements
     '''
    return [[bi.index,bi.index,0,0,0,-np.sign(bi.spin)*dS] for bi in basis]


def GrahamSchmidt(a,b):
    '''
    Simple orthogonalization of two vectors, returns orthonormalized vector
    args: a,b -- np.array of same length
    returns: tmp -- numpy array of same size, orthonormalized to the b vector
    '''
    tmp = a - np.dot(a,b)/np.dot(b,b)*b
    return tmp/np.linalg.norm(tmp)


def region(num):
    '''
    Generate a symmetric grid of points in number of lattice vectors. The tacit assumption is a 3 dimensional lattice
    args: num -- integer--grid will have size 2*num+1 in each direction
    returns numpy array of size ((2*num+1)**3,3) with centre value of first entry of (-num,-num,-num),...,(0,0,0),...,(num,num,num)
    '''
    num_symm = 2*num+1
    return np.array([[int(i/num_symm**2)-num,int(i/num_symm)%num_symm-num,i%num_symm-num] for i in range((num_symm)**3)])


        


# 
    

