#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 08:37:11 2020

@author: Rodrigo
"""

#%%
import matplotlib.pyplot as plt 

from parameterSettings import geometry_settings

from functions.manageDicom import readDicom
from functions.dataPreProcess import dataPreProcess
from functions.initialConfig import initialConfig

from functions.FBP import FDK as FBP


#%% Call function for initial configurations

libFiles = initialConfig()
    
#%% Create a DBT geometry
    
geo = geometry_settings()
geo.GE()

#%% Get DICOM data

dcmPath = "/home/rodrigo/Downloads/imgs/"
proj = readDicom(dcmPath,geo)

proj = dataPreProcess(proj, geo)

#%% Set specific recon parameters

filterType = 'FBP';          # Filter type: 'BP', 'FBP'
cutoff = 0.75;               # Percentage until cut off frequency (FBP)

#%%

print("Starting reconstruction...")

vol = FBP(proj, geo, filterType, cutoff, libFiles)


plt.figure()
plt.title('Reconstructed slice')
plt.imshow(vol[:,:,50] , cmap=plt.get_cmap('gist_gray'))

# print("Starting projection...")
# proj = projectionDDb(vol, geo, libFiles) 
# plt.figure()
# plt.title('Projected volume')
# plt.imshow(proj[:,:,4] , cmap=plt.get_cmap('gist_gray'))


