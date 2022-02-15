#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 18:46:40 2020

@author: Rodrigo
"""

import numpy as np
import numpy.ctypeslib as ctl

from .utilities import findAndLoadLibray
from .utilities import geoAsNp

def backprojectionDDb(proj, geo, libFiles):
    
    
    # Check if the input is in proper size
    if not (proj.shape[0] == geo.nv and proj.shape[1] == geo.nu and proj.shape[2] == geo.nProj):
        raise ValueError('First argument needs to have the same number of rows, cols and slices as in the configuration file.')
    
    
    # Find and library    
    lib = findAndLoadLibray(libFiles, 'backprojectionDDb')
    
    
    backprojectionDDb_lib = lib.backprojectionDDb_lib
    
    backprojectionDDb_lib.argtypes = [ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                               ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                               ctl.ndpointer(np.float32, flags='aligned, c_contiguous')]
    
    
    # Transform geo class in numpy array
    geoNp = geoAsNp(geo) 
    
    
    vol_transp = np.empty([geo.nz, geo.nx, geo.ny], dtype=np.float64)
    
    proj_transp = np.transpose(proj, (2, 1, 0)).copy()
    
    backprojectionDDb_lib (proj_transp, vol_transp, geoNp)
    
    vol = np.transpose(vol_transp, (2, 1, 0)).copy()
    
    
    return vol





















