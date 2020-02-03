#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 07:43:53 2020

@author: Rodrigo
"""

import numpy as np
from functions.backprojectionDDb import backprojectionDDb

def FDK(proj, geo, filterType, cutoff, libFiles):
    
    if filterType == 'FBP':
        print("----------------\nStarting FBP... \n\n")
        proj = filterProj(proj, geo, cutoff)
    elif filterType == 'BP':
        print("----------------\nStarting BP... \n\n")
    else:
        raise ValueError('Unknown filter type.')
        
    vol = backprojectionDDb(proj, geo, libFiles)
    
    return vol
    
    

def filterProj(proj, geo, cutoff):


    filteredProj = np.empty(proj.shape)
    
    us = np.linspace(geo.nu-1, 0, geo.nu) * geo.du
    vs = np.linspace(-(geo.nv-1)/2, (geo.nv-1)/2, geo.nv) * geo.dv
        
    # Detector Coordinate sytem in (mm)
    uCoord, vCoord = np.meshgrid(us, vs)
    
    # Compute weighted projections (Fessler Book Eq. (3.10.6))
    weightFunction = geo.DSO / np.sqrt((geo.DSD**2)+(vCoord**2)+(uCoord**2))
                       
    # Apply weighting function on each proj
    for i in range(geo.nProj):
        filteredProj[:,:,i] = proj[:,:,i] * weightFunction
    
    
    # Increase filter length to two times nv to decrease freq step
    h_Length = int(2**np.ceil(np.log2(np.abs(2*geo.nv))))
    
    # Builds ramp filter in space domain
    ramp_kernel = ramp_builder(h_Length)
    
    # Window filter in freq domain
    H_filter = filter_window(ramp_kernel, h_Length, cutoff)
    
    # Replicate to all colluns to build a 2D filter kernel
    H_filter = np.transpose(np.tile(H_filter, (geo.nu,1)))
    
    # Proj in freq domain
    H_proj = np.zeros([h_Length,geo.nu])
    
    # Filter each projection
    for i in range(geo.nProj):
         
       H_proj[0:geo.nv,:] = filteredProj[:,:,i]
    
       # Fourier transfor in projections
       H_proj = np.fft.fftshift(np.fft.fft(H_proj, axis=0))  
        
       # Multiplication in frequency = convolution in space
       H_proj = H_proj * H_filter
        
       # Inverse Fourier transfor
       H_proj = np.real(np.fft.ifft(np.fft.ifftshift(H_proj), axis=0))
    
       filteredProj[:,:,i] = H_proj[0:geo.nv,:]
    
    return filteredProj
    
    
    

## Function Ramp Filter
"""
The function builds Ramp Filter in space domain
Reference: Jiang Hsieh's book (second edition,page 73) Eq. 3.29
Reference: Fessler Book Eq.(3.4.14)
"""
def ramp_builder(h_Length):

    n = np.linspace(-h_Length/2, (h_Length/2)-1, h_Length)
    h = np.zeros(n.shape)
    h[int(h_Length/2)] = 1/4            # Eq. 3.29
    odd = np.mod(n,2) == 1              # Eq. 3.29
    h[odd] = -1 / (np.pi * n[odd])**2   # Eq. 3.29
    
    return h


## Function Window Ramp Filter
"""
The function builds Ramp Filter apodizided in frequency domain
Reference: Fessler Book and MIRT
"""
def filter_window(ramp_kernel, h_Length, cutoff):

    # Bring filter to freq domain
    H_ramp = np.abs(np.fft.fftshift(np.fft.fft(ramp_kernel))) 
    
    w = np.round(h_Length * cutoff) # Cut off freq
    n = np.linspace(-h_Length/2, (h_Length/2)-1, h_Length)
    H_window = 0.5 * (1 + np.cos(2 * np.pi * n /w)) # Hanning filter
    H_window = H_window * (np.abs(n) < w / 2) # Apply cut off freq
    
    H_filter = H_ramp * H_window # Apply window filter
    
    return H_filter
