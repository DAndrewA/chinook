#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 19:26:02 2017

@author: rday

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
import ubc_tbarpes.wigner as Wlib
import ubc_tbarpes.orbital as olib
import ubc_tbarpes.operator_library as oplib


def SK_coeff(o1,o2,R12,Vcoeff,renorm,offset,tol):
    if Vcoeff is not None:
        tmp = str(o1.atom)+str(o2.atom)+str(o1.label[0])+str(o2.label[0])+str(o1.label[1])+str(o2.label[1])
        tmp2 = str(o2.atom)+str(o1.atom)+str(o2.label[0])+str(o1.label[0])+str(o2.label[1])+str(o1.label[1])
#        print(tmp,tmp2)
        eta = Vcoeff
        
    E = {}
    SK_co = 0.0
    rad = np.linalg.norm(R12) #switched the two around 04/19/16 (o1 - o2 was old version)
    if o1.index == o2.index and rad<0.001:
        try:
            tmp=str(o1.atom)+str(o1.label) #in some cases, there may be different on-site energy for different orbitals of same l. So try this first, if not, then just do the l energy#format of orbital string must be ILO -- I = index, L = angular momentum (s,p,d,..), O = orientation (x,xy,ZR,...)
            SK_co = eta[tmp]-offset
        except KeyError:
            tmp = str(o1.atom)+o1.label[:2]
            SK_co = eta[tmp] -offset
    if rad>0.001:   
        l12 = (str(o1.l)+str(o2.l))
        l,m,n= np.array(R12)/rad  #switched o1-o2 to o2-o1 04/20/16
        if l12=="00":
#                    osS = o1.orbital[:3]+o2.orbital[:3]+"S"
            if tmp+"S" in eta:
                
                V = eta[tmp+"S"]
#                        V = eta[tmp+"S"]
            elif tmp2+"S" in eta:
                V = eta[tmp2+"S"]
#                        V = eta[tmp2+"S"]
            else:
                V = 0.0
            E["0-0"]= V
    
        elif l12=="01" or l12=="10":  

#                    osS = o1.orbital[:3]+o2.orbital[:3]+"S"
#                if eta[tmp+"S"]:
            if tmp+"S" in eta:
                V = eta[tmp+"S"]
#                        V = eta[tmp+"S"]
            elif tmp2+"S" in eta:
                V = -1*eta[tmp2+"S"] #for Vps, need to use -1*Vsp
            else:
                V = 0.0
            E["0-1x"],E["1x-0"] = [l*V]*2
            E["0-1y"],E["1y-0"] = [m*V]*2
            E["0-1z"],E["1z-0"] = [n*V]*2

         
        elif l12=="11":
            if tmp+"S" in eta:
                V = [eta[tmp+"S"],eta[tmp+"P"]]
#                        V = eta[tmp+"S"],eta[tmp+"P"]
            elif tmp2+"S" in eta:
                V = [eta[tmp2+"S"],eta[tmp2+"P"]]
#                        V = eta[tmp2+"S"],eta[tmp2+"P"]   
            else:
                V = [0.0,0.0]

            E["1x-1x"] = l**2*V[0] + (1.0-l**2)*V[1]
            E["1y-1y"] = m**2*V[0] + (1.0-m**2)*V[1]
            E["1z-1z"] = n**2*V[0] + (1.0-n**2)*V[1]
            E["1x-1y"],E["1y-1x"] = [l*m*(V[0] - V[1])]*2
            E["1x-1z"],E["1z-1x"] = [l*n*(V[0] - V[1])]*2
            E["1y-1z"],E["1z-1y"] = [m*n*(V[0] - V[1])]*2 
            
        elif l12 == "02" or l12 == "20":
#                    osS = o1.orbital[:3]+o2.orbital[:3]+"S"
            if tmp+"S" in eta:
                V = eta[tmp+"S"]
            elif tmp2+"S" in eta:
                V = eta[tmp2+"S"]
            else:
                V = 0.0
                    
            E["0-2xy"],E["2xy-0"] = [np.sqrt(3)*l*m*V]*2
            E["0-2yz"],E["2yz-0"] = [np.sqrt(3)*m*n*V]*2
            E["0-2xz"],E["2xz-0"] = [np.sqrt(3)*l*n*V]*2
            E["0-2XY"],E["2XY-0"] = [np.sqrt(3)/2.0*(l**2-m**2)*V]*2
            E["0-2ZR"],E["2ZR-0"] = [(n**2-(l**2+m**2)/2.0)*V]*2
                                
        elif l12== ("12") or l12 == "21":
#                    osS = o1.orbital[:3]+o2.orbital[:3]+"S"
#                    osP = o1.orbital[:3]+o2.orbital[:3]+"P"
                
            if tmp+"S" in eta:
                V = [eta[tmp+"S"],eta[tmp+"P"]]
            elif tmp2+"S" in eta:
                V = -1*[eta[tmp2+"S"],eta[tmp2+"P"]]
            else:
                V = [0.0,0.0]
                    
            E["1x-2xy"],E["2xy-1x"] = [np.sqrt(3)*l**2*m*V[0]+m*(1-2*l**2)*V[1]]*2
            E["1x-2yz"],E["2yz-1x"] = [np.sqrt(3)*l*m*n*V[0]-2*l*m*n*V[1]]*2
            E["1x-2xz"],E["2xz-1x"] = [np.sqrt(3)*l**2*n*V[0]+n*(1-2*l**2)*V[1]]*2
            
            E["1y-2xy"],E["2xy-1y"] = [np.sqrt(3)*m**2*l*V[0]+l*(1-2*m**2)*V[1]]*2
            E["1y-2xz"],E["2xz-1y"] = [np.sqrt(3)*l*m*n*V[0]-2*l*m*n*V[1]]*2
            E["1y-2yz"],E["2yz-1y"] = [np.sqrt(3)*m**2*n*V[0]+n*(1-2*m**2)*V[1]]*2                    
            
            E["1z-2yz"],E["2yz-1z"] = [np.sqrt(3)*n**2*m*V[0]+m*(1-2*n**2)*V[1]]*2
            E["1z-2xy"],E["2xy-1z"] = [np.sqrt(3)*l*m*n*V[0]-2*l*m*n*V[1]]*2
            E["1z-2xz"],E["2xz-1z"] = [np.sqrt(3)*n**2*l*V[0]+l*(1-2*n**2)*V[1]]*2
                
            E["1x-2XY"],E["2XY-1x"] = [np.sqrt(3)/2*l*(l**2-m**2)*V[0]-l*(1-l**2+m**2)*V[1]]*2
            E["1y-2XY"],E["2XY-1y"] = [np.sqrt(3)/2*m*(l**2-m**2)*V[0]-m*(1-m**2+l**2)*V[1]]*2
            E["1z-2XY"],E["2XY-1z"] = [np.sqrt(3)/2*n*(l**2-m**2)*V[0] - n*(l**2-m**2)*V[1]]*2
                
            E["1x-2ZR"],E["2ZR-1x"] = [l*(n**2-(l**2+m**2)/2.0)*V[0]-np.sqrt(3)*l*n**2*V[1]]*2
            E["1y-2ZR"],E["2ZR-1y"] = [m*(n**2-(l**2+m**2)/2.0)*V[0]-np.sqrt(3)*m*n**2*V[1]]*2
            E["1z-2ZR"],E["2ZR-1z"] = [n*(n**2-(l**2+m**2)/2.0)*V[0]+np.sqrt(3)*n*(l**2+m**2)*V[1]]*2
            
                
        elif l12== "22":

            if tmp+"S" in eta:
                V = [eta[tmp+"S"],eta[tmp+"P"],eta[tmp+"D"]]
            elif tmp2+"S" in eta:                       
                V = [eta[tmp2+"S"],eta[tmp2+"P"],eta[tmp2+"D"]]
            else:
                V = [0.0,0.0,0.0]     
            E["2xy-2xy"] = 3*l**2*m**2*V[0]+(l**2+m**2-4*l**2*m**2)*V[1]+(n**2+l**2*m**2)*V[2]
            E["2yz-2yz"] = 3*m**2*n**2*V[0]+(m**2+n**2-4*m**2*n**2)*V[1]+(l**2+m**2*n**2)*V[2]
            E["2xz-2xz"] = 3*n**2*l**2*V[0]+(n**2+l**2-4*n**2*l**2)*V[1]+(m**2+n**2*l**2)*V[2]

            E["2xy-2yz"],E["2yz-2xy"] = [3*l*m**2*n*V[0]+l*n*(1-4*m**2)*V[1]+l*n*(m**2-1)*V[2]]*2
            E["2xy-2xz"],E["2xz-2xy"] = [3*m*l**2*n*V[0]+m*n*(1-4*l**2)*V[1]+m*n*(l**2-1)*V[2]]*2
            E["2yz-2xz"],E["2xz-2yz"] = [3*l*n**2*m*V[0]+l*m*(1-4*n**2)*V[1]+l*m*(n**2-1)*V[2]]*2
                
            E["2xy-2XY"],E["2XY-2xy"] = [1.5*l*m*(l**2-m**2)*V[0]+2*l*m*(m**2-l**2)*V[1]+(l*m*(l**2-m**2)/2.0)*V[2]]*2
            E["2xy-2ZR"],E["2ZR-2xy"] = [np.sqrt(3)*((l*m*(n**2-(l**2+m**2)/2.0))*V[0]-2*l*m*n**2*V[1]+(l*m*(1+n**2)/2)*V[2])]*2
            
            E["2yz-2XY"],E["2XY-2yz"] = [1.5*m*n*(l**2-m**2)*V[0]-m*n*(1+2*(l**2-m**2))*V[1]+(m*n*(1+(l**2-m**2)/2.0))*V[2]]*2
            E["2yz-2ZR"],E["2ZR-2yz"] = [np.sqrt(3)*((m*n*(n**2-(l**2+m**2)/2.0))*V[0]+m*n*(l**2+m**2-n**2)*V[1]-(m*n*(l**2+m**2)/2)*V[2])]*2     
            
            E["2xz-2XY"],E["2XY-2xz"] = [1.5*n*l*(l**2-m**2)*V[0]+n*l*(1-2*(l**2-m**2))*V[1]-(n*l*(1-(l**2-m**2)/2.0))*V[2]]*2
            E["2xz-2ZR"],E["2ZR-2xz"] = [np.sqrt(3)*((l*n*(n**2-(l**2+m**2)/2.0))*V[0]+l*n*(l**2+m**2-n**2)*V[1]-(l*n*(l**2+m**2)/2)*V[2])]*2
                
            E["2XY-2XY"] = 0.75*(l**2-m**2)**2*V[0]+(l**2+m**2-(l**2-m**2)**2)*V[1]+(n**2+(l**2-m**2)**2/4.0)*V[2]
            E["2ZR-2ZR"] = (n**2-(l**2+m**2)/2.0)**2*V[0]+3*n**2*(l**2+m**2)*V[1]+0.75*(l**2+m**2)**2*V[2]
            E["2XY-2ZR"],E["2ZR-2XY"] = [np.sqrt(3)*((l**2-m**2)*(n**2-(l**2+m**2)/2.0)*V[0]/2.0+n**2*(m**2-l**2)*V[1]+((1+n**2)*(l**2-m**2)/4.0)*V[2])]*2

           ####DISTORTION NOT WORKING CORRECTLY YET!
#           if type(o1.Dmat)==np.ndarray or type(o2.Dmat)==np.ndarray:
#                ldict = {0:{'0':0},1:{'1x':0,'1y':1,'1z':2},2:{'2xz':0,'2yz':1,'2xy':2,'2ZR':3,'2XY':4}}
#                cdict = {}
#                L1 = np.identity(2*o1.l+1)
#                L2 = np.identity(2*o2.l+1)
#                if o1.Dmat is not None:
#                    L1 = o1.Dmat
#                if o2.Dmat is not None:
#                    L2 = o2.Dmat
#                tmat = np.zeros((2*(o1.l+o2.l+1),2*(o1.l+o2.l+1)),dtype=complex)
#                tmat[:2*o1.l+1,:2*o1.l+1] = L1
#                tmat[2*o1.l+1:,2*o1.l+1:] = L2
#                H = np.zeros(np.shape(tmat))
#                for e in E:
#                    st = e.split('-')
#                    cdict[ldict[o1.l][st[0]]] = st[0]
#                    cdict[ldict[o2.l][st[1]]+o1.l*2+1] = st[1]
#                    H[ldict[o1.l][st[0]],ldict[o2.l][st[1]]+o1.l*2+1] = E[e]
#                E = {}
#                Ht = np.dot(np.linalg.inv(tmat),np.dot(H,tmat))
#                for i in range(np.shape(Ht)[0]):
#                    for j in range(np.shape(Ht)[1]):
#                        if abs(Ht[i,j])>10**-10:
#                            print(i,j,Ht[i,j])
#                        E[cdict[i]+'-'+cdict[j]] = Ht[i,j]
#                

     
        SK_co = E[o1.label[1:]+'-'+o2.label[1:]]
    SK_co/=renorm
    if abs(SK_co)<tol:
        SK_co=0.0
    return SK_co#,E


def eta_SK(eta,d,l,lp,rd):
    hb2m = 7.62
    if lp==2:
        if l==2:
            return eta*hb2m*rd**3/d**5
        else:
            return eta*hb2m*rd**1.5/d**3.5
    else:
        return eta*hb2m/d**2
    
def eta_dict_to_SK(eta_dict,rd,d):
    VSK = {}
    for di in eta_dict:
        try:
            l = int(di[4])
            lp = int(di[5])
            VSK[di] = eta_SK(eta_dict[di],d,l,lp,rd)
        except IndexError:
            VSK[di] = eta_dict[di]
    return VSK


def Vmat(l1,l2,V):
    '''
    For Slater-Koster matrix element generation, a potential matrix is sandwiched in between the
    two bond-rotating Dmatrices. It should be of the shape 2*l1+1 x 2*l2+1, and have the Vl,l',D
    terms along the 'diagonal'-- a concept that is only well defined for a square matrix. For mismatched
    angular momentum channels, this turns into a diagonal square matrix of dimension min(2*l1+1,2*l2+1) 
    centred  along the larger axis. For channels where the orbital angular momentum change involves a change
    in parity, the potential should change sign, as per Slater Koster's original definition from 1954. 
    args:
        l1,l2: int orbital angular momentum of initial and final states
        V: numpy array of float -- length should be min(l1,l2)*2+1
    return:
        Vm: numpy array of float, shape 2*l1+1 x 2*l2+1
    '''
    if l1==0 or l2==0:
        Vm = np.zeros(max(2*l1+1,2*l2+1))
        Vm[int(len(Vm)/2)] = V[0]
    else:
        Vm = np.zeros((2*l1+1,2*l2+1))
        lmin = min(l1,l2)
        lmax = max(l1,l2)
        Vvals = np.identity(2*lmin+1)*np.array([V[abs(i-lmin)] for i in range(2*lmin+1)])
        if l2>l1:
            Vm[:,lmax-lmin:lmax-lmin+2*lmin+1] = Vvals
        else:
            Vm[lmax-lmin:lmax-lmin+2*lmin+1,:] = Vvals
    if l1>l2 and (np.mod(l2-l1,2)!=0):
        Vm*=-1
    return Vm

def SK_Wig_Y(l1,l2):
    '''
    Generate a matrix function that can be applied to any hopping-combination of the orbitals labelled by
    l1,l2. The output, a lambda function can then be addressed with a given bond and potential, generating
    the hoppings. Right now, this is still in the basis of spherical harmonics. Need to prepend by the local
    projection array on either side. Then to implement in generating a Hamiltonian, one would generate the relevant
    lambda functions for the different angular momentum channels involved, then evaluate the matrix for a specific
    choice of hopping. In principal, this means we could just iterate over atomic sites rather than all combinations
    of orbitals, and then just pull all the matrix elements from the arrays evaluated at these points and fill them
    in. This could be very efficient method, and more importantly, generic for all angular momentum channels!
    To be completed 26/9/2018
    '''
    

    return lambda A,B,V: np.dot(np.conj(Wlib.WignerD(l1,A,B,0).T),np.atleast_2d(np.dot(Vmat(l1,l2,V),Wlib.WignerD(l2,A,B,0))))


def SK_cub(Ymats,l1,l2):
    def SK_build(A,B,V):
    
        return np.dot(np.dot(np.conj(Ymats[0].T),np.conj(Wlib.WignerD(l1,A,B,0).T)),np.atleast_2d(np.dot(Vmat(l1,l2,V),np.dot(Wlib.WignerD(l2,A,B,0),Ymats[1]))))
    return lambda A,B,V:SK_build(A,B,V)

def SK_full(basis):
    
    '''
    Generate a dictionary of lambda functions which take as keys the atom,orbital for both first and second element 
    '''
    SK_funcs = {}
    Ymats = oplib.Yproj(basis)
    for yi in Ymats:
        for yj in Ymats:
            Y = [Ymats[yi],Ymats[yj]]
            SK_funcs[(yi[0],yi[1],yj[0],yj[1])] = SK_cub(Y,yi[1],yj[1])
    return SK_funcs


        
    



    


#if __name__=="__main__":
#    Vo=1./6.25
#    off = 0.5#1.1722
#    scale = np.sqrt(2)*5.48
#    dxy = -0.0
#    
#    SK = {"052xy":-off+dxy,"152xy":-off+dxy,"052xz":-off,"152xz":-off,"052yz":-off,"152yz":-off,"052XY":-off-3,"152XY":-off-3,"052ZR":-off-3,"152ZR":-off-3,"015522S":1.5*Vo,"015522P":-1.0*Vo,"015522D":0.25*Vo,"115522S":1.5*Vo/scale,"115522P":-1.0*Vo/scale,"115522D":0.25*Vo/scale,"005522S":1.5*Vo/scale,"005522P":-1.0*Vo/scale,"005522D":0.25*Vo/scale}
#    pos1 = np.zeros(3)
#    pos2 = np.array([3.875,0,0])
#    xz1 = olib.orbital(0,0,"52xz",pos1,77)#,orient=[np.array([0,0,1]),90*np.pi/180.])
#    yz1 = olib.orbital(0,1,"52yz",pos1,77)#,orient=[np.array([0,0,1]),90*np.pi/180.])
#    xy1 = olib.orbital(0,2,"52xy",pos1,77)#,orient=[np.array([0,0,1]),90*np.pi/180.])
#    xz2 = olib.orbital(1,3,"52xz",pos2,77)#,orient=[np.array([0,0,1]),-90*np.pi/180.])
#    yz2 = olib.orbital(1,4,"52yz",pos2,77)#,orient=[np.array([0,0,1]),-90*np.pi/180.])
#    xy2 = olib.orbital(1,5,"52xy",pos2,77)#,orient=[np.array([0,0,1]),-90*np.pi/180.])
    
#    V,E = SK_coeff(xz1,yz2,SK)
#    os = [xz1,yz1,xy1,xz2,yz2,xy2]
#    for o in os:
#        for o2 in os:
#            V,_ = SK_coeff(os[0],os[4],SK)
#            if abs(V)!=0:
#                print(o.label,o2.label,V)
#    
#    