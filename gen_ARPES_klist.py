# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 13:57:58 2022

@author: Andrew

Function to generate a case.klist_band file that encompasses a small region 
around a selected point in the Brillouen Zone, for use in a bands calculation
by WIEN2K for the purpose of being used with the Chinook library.
"""

import numpy as np

def gen_ARPES_klist(case,center,extent,subdivisions):
    '''
    This function will generate a case.klist_band file for the use in WIEN2K that
    follows the format required by datacube in the Chinook library's ARPES experiment.
    '''
    print('bah')

