#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:34:24 2019

@author: rodrigo
"""

import numpy as np
from scipy.ndimage.filters import uniform_filter1d
from skimage.filters import threshold_otsu

def dataPreProcess(proj, geo, flagtransfIntensity=True, flagCropProj=True):
    
    if flagCropProj:
        proj = cropProj(proj)
    
    if flagtransfIntensity:
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
    indX = np.argmin(secondDv)  
    
    # Alternative method to compare
    # Otsu
    thresh = threshold_otsu(proj[:,:,proj.shape[2]//2 - 1])
    
    # Get mask
    mask = proj[:,:,proj.shape[2]//2 - 1] < thresh
    
    # Weight bounds
    mask_w = np.sum(mask, 0) > 0
    res = np.where(mask_w == True)
    w_min, w_max = res[0][0], res[0][-1]
    
    # If the difference is big, take the min (most left)
    if np.abs(indX - w_min) > 200:
        indX = np.min([indX, w_min])
    
    indX -= Gap
        
    proj_crop = proj[:,indX::,:]
    		    
    return proj_crop


def transfIntensity(proj):
    
     # Transform intensity image in attenuation coefficients
    proj = -np.log(proj/np.max(proj)); 
    
    return proj
    
