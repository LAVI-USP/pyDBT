#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:29:14 2022

@author: rodrigo
"""

import numpy as np
import numpy.ctypeslib as ctl
import ctypes
import warnings

from .utilities import findAndLoadLibray
from .utilities import geoAsNp

##############################################################################
#                         BackProjection                                     #
##############################################################################

def backprojectionDD(proj, geo, proj_num, libFiles):
    
    operator_type = 'backprojectionDD'
    
    vol = _backprojection(proj, geo, proj_num, libFiles, operator_type)
    
    return vol

def backprojectionDDb(proj, geo, proj_num, libFiles):
    
    operator_type = 'backprojectionDDb'
    
    vol = _backprojection(proj, geo, proj_num, libFiles, operator_type)
        
    return vol

def backprojectionDDb_cuda(proj, geo, proj_num, libFiles):
    
    operator_type = 'backprojectionDDb_cuda'
    
    vol = _backprojection(proj, geo, proj_num, libFiles, operator_type)
    
    return vol
    
def _backprojection(proj, geo, proj_num, libFiles, operator_type):
    
    # Check if the input is in proper size
    if not (proj.shape[0] == geo.nv and proj.shape[1] == geo.nu and proj.shape[2] == geo.nProj):
        raise ValueError('First argument needs to have the same number of rows, cols and slices as in the configuration file.')
    
    check_parameters(operator_type, proj_num, geo)
    
    # Find and library    
    lib = findAndLoadLibray(libFiles, operator_type)
    
    
    backprojection = getattr(lib, operator_type + '_lib')
    
    backprojection.argtypes = [ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                               ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                               ctl.ndpointer(np.float32, flags='aligned, c_contiguous'),
                               ctypes.c_long]
    
    
    # Transform geo class in numpy array
    geoNp = geoAsNp(geo) 
    
    
    vol_transp = np.empty([geo.nz, geo.nx, geo.ny], dtype=np.float64)
    
    proj_transp = np.transpose(proj, (2, 1, 0)).copy()
    
    backprojection(proj_transp, vol_transp, geoNp, proj_num)
    
    vol = np.transpose(vol_transp, (2, 1, 0)).copy()
    
    
    return vol
    

##############################################################################
#                           Projection                                       #
##############################################################################

def projectionDD(vol, geo, proj_num, libFiles):
    
    operator_type = 'projectionDD'
    
    proj = _projection(vol, geo, proj_num, libFiles, operator_type)
    
    return proj

def projectionDDb(vol, geo, proj_num, libFiles):
    
    operator_type = 'projectionDDb'
    
    proj = _projection(vol, geo, proj_num, libFiles, operator_type)
    
    return proj

def projectionDDb_cuda(vol, geo, proj_num, libFiles):
    
    operator_type = 'projectionDDb_cuda'
    
    proj = _projection(vol, geo, proj_num, libFiles, operator_type)
    
    return proj

def _projection(vol, geo, proj_num, libFiles, operator_type):
    
    # Check if the input is in proper size
    if not (vol.shape[0] == geo.ny and vol.shape[1] == geo.nx and vol.shape[2] == geo.nz):
        raise ValueError('First argument needs to have the same number of rows, cols and slices as in the configuration file.')
       
    check_parameters(operator_type, proj_num, geo)
        
    # Find and library    
    lib = findAndLoadLibray(libFiles, operator_type)
    
    projection = getattr(lib, operator_type + '_lib')
    
    projection.argtypes = [ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                           ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                           ctl.ndpointer(np.float32, flags='aligned, c_contiguous'),
                           ctypes.c_long]
    
    
    # Transform geo class in numpy array
    geoNp = geoAsNp(geo)
    
    
    proj_transp = np.empty([geo.nProj, geo.nu, geo.nv], dtype=np.float64)
    
    vol_transp = np.transpose(vol, (2, 1, 0)).copy()
        
    projection(vol_transp, proj_transp, geoNp, proj_num)
    
    proj = np.transpose(proj_transp, (2, 1, 0)).copy()
    
    return proj

##############################################################################
#                           Auxiliar functions                               #
##############################################################################


def check_parameters(operator_type, proj_num, geo):
    
    if proj_num < -1 and proj_num >= geo.nProj:
        raise ValueError('Projection number needs to be between 0-(geo.nProj-1). If you want to operate on all projections, set it to -1')
        
    if geo.detAngle != 0 and 'DDb' in operator_type:
        warnings.warn('Your detector rotates, its not recomended to used branchless projectors. Use DD only')
        
    if (geo.x_offset != 0 or geo.y_offset != 0) and ('DDb' in operator_type):
        raise ValueError('Branchless operators dont work with volume offsets')
        
    return 
    
    
    
    
