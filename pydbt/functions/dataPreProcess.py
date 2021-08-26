#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:34:24 2019

@author: rodrigo
"""

import numpy as np
from scipy.ndimage.filters import uniform_filter1d

def dataPreProcess(proj, geo, flagCropProj=True):
    
    if flagCropProj:
        proj = cropProj(proj)
    
    proj = transfIntensity(proj)
    
    # Modifies parameters based on segmentation
    geo.nu = proj.shape[1]  # Number of pixels (columns)
    
    return proj
    


def cropProj(proj):
    
    Gap = 20;
    
    # Horizontal Profile
    vertProj = np.sum(proj[:,:,proj.shape[2]//2 - 1], axis=0)	
    
    # Smooth the signal and take the first derivative
    firstDv = np.gradient(uniform_filter1d(vertProj, size=100))
    
    # Smooth the signal and take the second derivative
    secondDv = np.gradient(uniform_filter1d(firstDv, size=100))
    
    # Takes its min second derivative
    indX = np.argmin(secondDv) - Gap 
        
    proj_crop = proj[:,indX::,:]
    		    
    return proj_crop


def transfIntensity(proj):
    
     # Transform intensity image in attenuation coefficients
    proj = -np.log(proj/np.max(proj)); 
    
    return proj
    
