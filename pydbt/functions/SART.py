#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 11:29:43 2022

@author: rodrigo
"""

import numpy as np
import time

from .projection_operators import backprojectionDDb, projectionDDb

def SART(proj, geo, nIter, libFiles):
       
    highestValue = (2**16) - 1 
    
    # Initial estimated data 
    reconData3d = np.zeros([geo.ny, geo.nx, geo.nz]) 
    
    # Pre calculation of Projection normalization
    proj_norm = projectionDDb(np.ones([geo.ny, geo.nx, geo.nz]), geo, -1, libFiles) 
    proj_norm[proj_norm == 0] = 1
    
    # Pre calculation of Backprojection normalization
    vol_norm = backprojectionDDb(np.ones([geo.nv, geo.nu, geo.nProj]), geo, -1, libFiles) 
    vol_norm[vol_norm == 0] = 1
    
    print('----------------\nStarting SART Iterations... \n\n')
    
    # Start Iterations
    for n_iter in range(nIter[-1]):
        
        startTime = time.time()
        
        for proj_num in range(geo.nProj):
        
            # Error between raw data and projection of estimated data  
            proj_diff = proj - projectionDDb(reconData3d, geo, proj_num, libFiles)  
            
            # Projection normalization
            proj_diff = proj_diff / proj_norm  
            proj_diff[np.isnan(proj_diff)] = 0 
            proj_diff[np.isinf(proj_diff)] = 0 
            
            upt_term = backprojectionDDb(proj_diff, geo, proj_num, libFiles) 
            upt_term = upt_term / vol_norm  # Volume normalization
            upt_term[np.isnan(upt_term)] = 0 
            upt_term[np.isinf(upt_term)] = 0 
        
            # Updates the previous estimation
            reconData3d = reconData3d + upt_term   
            
            endTime = time.time()
                
            # Truncate to highest and minimum value        
            reconData3d[reconData3d > highestValue] = 0 
            reconData3d[reconData3d < 0] = 0 
        
        print('Itaration %d Time: %.2fs\n\n' % (n_iter , endTime-startTime)) 
    
    return reconData3d