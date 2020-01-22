#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 08:39:59 2019

@author: rodrigo
"""

import pydicom
import numpy as np
from os import walk

def readDicom(path, geo):
    
    dcmFiles = []
    for (_, _, filenames) in walk(path):
        for file in filenames:
            if file.endswith(".dcm"):
                dcmFiles.append(file)
        break
    
    # Test if list is empty
    if not dcmFiles:    
        raise ValueError('No DICOM files found in the specified path.')
    
    proj = np.zeros([geo.nv,geo.nu,geo.nProj],dtype='uint16')
    
    
    for f in dcmFiles:
        nProj = int(f.split('.')[0])
        proj_dcmH = pydicom.dcmread(path + f)
        proj[:,:,nProj]  = proj_dcmH.pixel_array
    
    return proj