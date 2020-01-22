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
from functions.projectionDDb import projectionDDb
from functions.backprojectionDDb import backprojectionDDb
from functions.initialConfig import initialConfig

#%% Call function for initial configurations

libFiles = initialConfig()
    
#%% Create a DBT geometry
    
geo = geometry_settings()
geo.GE()

#%% Get DICOM data

dcmPath = "/home/rodrigo/Downloads/imgs/"
proj = readDicom(dcmPath,geo)

proj = dataPreProcess(proj, geo)

#%%

print("Starting reconstruction...")

vol = backprojectionDDb(proj, geo, libFiles)
plt.figure()
plt.title('Reconstructed slice')
plt.imshow(vol[:,:,50] , cmap=plt.get_cmap('gist_gray'))

print("Starting projection...")

proj = projectionDDb(vol, geo, libFiles) 
plt.figure()
plt.title('Projected volume')
plt.imshow(proj[:,:,4] , cmap=plt.get_cmap('gist_gray'))


