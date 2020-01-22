#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:34:24 2019

@author: rodrigo
"""

import numpy as np

def dataPreProcess(proj, geo):
    
    proj = cropProj(proj)
    
    proj = transfIntensity(proj)
    
    # Modifies parameters based on segmentation
    geo.nu = proj.shape[1]  # Number of pixels (columns)
    
    return proj
    


def cropProj(proj):
    
    Gap = 20;
    
    maxValue = np.max(proj)
    
    vertProj = np.sum(maxValue-proj[:,:,proj.shape[2]//2 - 1], axis=0)	# Horizontal Profile
    
    Ind = np.argmax(np.diff(vertProj[10:])) - Gap   # Smooth the signal and takes its max positive derivative
    
    proj_crop = proj[:,Ind:,:]
    		    
    return proj_crop


def transfIntensity(proj):
    
     # Transform intensity image in attenuation coefficients
    proj = -np.log(proj/np.max(proj)); 
    
    return proj
    
