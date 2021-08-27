#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:35:00 2020

@author: rodrigo
"""
 
import pathlib

def initialConfig():

    workDir = pathlib.Path().absolute()
    
    outDir = workDir / 'output'
    
    if not(outDir.exists() and outDir.is_dir()):
        pathlib.Path(outDir).mkdir()
                
    libFiles = findLibraries(workDir)
    
    return libFiles
        

def findLibraries(workDir):
    
    buildDir = workDir / '../build' 
    
    if not(buildDir.exists() and buildDir.is_dir()):
        
        buildDir = workDir / '../../build' 
        
        if not(buildDir.exists() and buildDir.is_dir()):
        
            raise ValueError('Cannot find the build folder. Make sure to run the setup.py or follow the instructions on the package Github.')
    
    libDir = []
    
    # Find libray directory
    for f in buildDir.iterdir():
        if f.is_dir():
            if 'lib.' in str(f):
                libDir = f
        
    # Test if libDir is empty
    if not libDir:    
        raise ValueError('Cannot find the lib folder. Make sure to run the setup.py or follow the instructions on the package Github.')
        
    libFiles = []
        
    # Find all .so files   
    for libFile in pathlib.Path(libDir).glob('*.so'):
        libFiles.append((str(libDir),str(libFile).split('/')[-1]))
        
        
    return libFiles

