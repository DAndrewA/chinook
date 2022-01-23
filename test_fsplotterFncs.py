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

import matplotlib as plt



case = 'data/FeTe_19nov21'


avec,a,b,c = fs.get_PLV_from_struct(case)
k_object = fs.create_kobject(case)

# load in QTL, comment out if not nessecary as very memory intensive
QTL,bands = fs.load_qtl(case)


# FIRST TEST:
    # are the lengths of the kpath and bands data the same?
if len(k_object.kcut) == len(bands[1]):
    print('Test 1: lengths are the same!')


#fig = plt.figure()


# plot bands from bands and k_object data:
