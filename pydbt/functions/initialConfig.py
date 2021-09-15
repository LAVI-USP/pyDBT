#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:35:00 2020

@author: rodrigo
"""
 
import pathlib

def initialConfig(buildDir=None, createOutFolder=True):

    # Build folder was not provided, try find it on up folders
    if buildDir == None:
        workDir = pathlib.Path().absolute()
        buildDir = findBuildDir(workDir)
    else:
        buildDir = pathlib.Path(buildDir)
    
    # Should I create output folder?
    if createOutFolder:
        outDir = workDir / 'output'
        
        if not(outDir.exists() and outDir.is_dir()):
            pathlib.Path(outDir).mkdir()
    
    # Find the compiled lib files
    libFiles = findLibraries(buildDir)
    
    return libFiles


def findBuildDir(workDir):
    
    buildDir = workDir / '../build' 
    
    if not(buildDir.exists() and buildDir.is_dir()):
        
        buildDir = workDir / '../../build' 
        
        if not(buildDir.exists() and buildDir.is_dir()):
        
            raise ValueError('Cannot find the build folder. Make sure to run the setup.py or follow the instructions on the package Github.')
    
    return buildDir
       

def findLibraries(buildDir):
    
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

