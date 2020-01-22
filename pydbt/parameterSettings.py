#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:19:22 2020

@author: Rodrigo
"""
     
class geometry_settings():
    
    # def __init__(self):
        
        # Intialize some empty variables for later on
        # self.nx = None
        # self.ny = None    
        # self.nz = None    
        # self.nu = None   
        # self.nv = None   
        # self.dx = None 
        # self.dy = None
        # self.dz = None
        # self.du = None   
        # self.dv = None
        # self.DSD = None                           
        # self.DSO = None                          
        # self.DDR = None                            
        # self.DSR = None          
        # self.DAG = None                
        # self.nProj = None  
        # self.tubeAngle = None
        # self.detAngle = None
        
    
    def GE(self):
         
        # Breast voxels density 
        self.nx = 1058    # number of voxels (columns)
        self.ny = 1978    # number of voxels (rows)
        self.nz = 107    # number of voxels (slices)
        
        # Detector panel pixel density
        self.nu = 2394    # number of pixels (columns)
        self.nv = 3062    # number of pixels (rows)
        
        # Single voxel real size (mm)
        self.dx = 0.1 
        self.dy = 0.1
        self.dz = 0.5
        
        # Single detector real size (mm)
        self.du = 0.1   
        self.dv = 0.1
        
        # X-ray source and detector distances
        self.DSD = 660                           # Distance from source to detector (mm)
        self.DSO = 580.5                         # Distance from source to the top of object (mm)
        self.DDR = 40                            # Distance from detector to pivot (mm)
        self.DSR = self.DSD - self.DDR           # Distance from source to pivot (mm)
        self.DAG = 22                            # Distance of Air Gap (mm)         
        
        
        # Number of Projections
        self.nProj = 9  
        
        # Angle settings (Degrees)
        self.tubeAngle = 25     # Tube Angle
        
        self.detAngle = 0       # Detector Angle
 