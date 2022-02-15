#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 13:36:50 2020

@author: Rodrigo
"""

import numpy as np
import numpy.ctypeslib as ctl

from .utilities import findAndLoadLibray
from .utilities import geoAsNp

def projectionDDb(vol, geo, libFiles):
    

    # Check if the input is in proper size
    if not (vol.shape[0] == geo.ny and vol.shape[1] == geo.nx and vol.shape[2] == geo.nz):
        raise ValueError('First argument needs to have the same number of rows, cols and slices as in the configuration file.')
        
        
    # Find and library    
    lib = findAndLoadLibray(libFiles, 'projectionDDb')
    
    
    projectionDDb_lib = lib.projectionDDb_lib
    
    projectionDDb_lib.argtypes = [ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                                  ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                                  ctl.ndpointer(np.float32, flags='aligned, c_contiguous')]
    
    
    # Transform geo class in numpy array
    geoNp = geoAsNp(geo)
    
    
    proj_transp = np.empty([geo.nProj, geo.nu, geo.nv], dtype=np.float64)
    
    vol_transp = np.transpose(vol, (2, 1, 0)).copy()
        
    projectionDDb_lib (vol_transp, proj_transp, geoNp)
    
    proj = np.transpose(proj_transp, (2, 1, 0)).copy()
    
    
    return proj

            
      
        

    