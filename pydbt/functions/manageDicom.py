#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 08:39:59 2019

@author: rodrigo
"""

import pydicom
import numpy as np
from pathlib import Path

def readDicom(path, geo):
    
    dcmFiles = [str(item) for item in Path(path).glob("*.dcm")]
    
    # Test if list is empty
    if not dcmFiles:    
        raise ValueError('No DICOM files found in the specified path.')
    
    proj = [None] * geo.nProj
    proj_header = [None] * geo.nProj
    
    for f in dcmFiles:
        nProj = int(f.split('/')[-1].split('.')[0])
        proj_header[nProj] = pydicom.dcmread(f)
        proj[nProj]  = proj_header[nProj].pixel_array
    
    proj = np.stack(proj, axis=-1).astype(np.uint16)
    
    return proj, proj_header