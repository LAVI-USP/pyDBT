#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:45:08 2020

@author: rodrigo
"""

import numpy.ctypeslib as ctl
import numpy as np


def findAndLoadLibray(libFiles, libname):
    
    for libFile in libFiles:
        if libname == libFile[1].split('.')[0]:
            return ctl.load_library(libFile[1], libFile[0])
        else:
            ValueError('Cannot find ' + libname + ' libray. Make sure it is compiled.' )

def geoAsNp(geo):
    
    geoList = []

    geoDict = geo.__dict__
    
    for key, value in geoDict.items():
        geoList.append(value)
    
    geoNp = np.asarray(geoList, dtype=np.float32)
    
    return geoNp
