#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 13:49:27 2020

@author: Rodrigo
"""

#%%
import matplotlib.pyplot as plt 
import numpy as np
import sys

sys.path.insert(1,'../')

from parameters.parameterSettings import geometry_settings

from functions.initialConfig import initialConfig
from functions.phantoms import phantom3d

from functions.projectionDD import projectionDD

from functions.FBP import FDK as FBP


#%% Call function for initial configurations

libFiles = initialConfig()
    
#%% Create a DBT geometry
    
geo = geometry_settings()
geo.SheppLogan()

#%% Gen Phantom

sheppLogan = phantom3d(n=geo.nx)
sheppLogan = sheppLogan[:,:,np.round(np.linspace(0,geo.nx-1,geo.nz)).astype(np.int32)]

#%% Project Phantom

print("Starting projection...")

proj = projectionDD(sheppLogan, geo, libFiles) 

plt.figure()
plt.title('Projected volume')
plt.imshow(proj[:,:,4] , cmap=plt.get_cmap('gist_gray'))

#%% Set specific recon parameters

nIter = [2];                 # Number of iterations (SIRT)
filterType = 'FBP';          # Filter type: 'BP', 'FBP'
cutoff = 0.75;               # Percentage until cut off frequency (FBP)

#%%

print("Starting reconstruction...")

#                       ## Uncomment to use ##
vol = FBP(proj, geo, filterType, cutoff, libFiles)

plt.figure()
plt.title('Reconstructed slice')
plt.imshow(vol[:,:,64] , cmap=plt.get_cmap('gist_gray'))
