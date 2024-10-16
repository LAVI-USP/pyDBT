#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:19:22 2020

@author: Rodrigo
"""
     
class geometry_settings():
    
    def Hologic(self):
         
        # Breast voxels density 
        self.nx = 1996    # number of voxels (columns)
        self.ny = 2457    # number of voxels (rows)
        self.nz = 78     # number of voxels (slices)
        
        # Detector panel pixel density
        self.nu = 1664    # number of pixels (columns)
        self.nv = 2048    # number of pixels (rows)
        
        # Single voxel real size (mm)
        self.dx = 0.112 
        self.dy = 0.112
        self.dz = 1
        
        # Single detector real size (mm)
        self.du = 0.14  
        self.dv = 0.14
        
        # X-ray source and detector distances
        self.DSD = 700                           # Distance from source to detector (mm)
        self.DSO = 597                           # Distance from source to the top of object (mm)
        self.DDR = 0                             # Distance from detector to pivot (mm)
        self.DSR = self.DSD - self.DDR           # Distance from source to pivot (mm)
        self.DAG = 25                            # Distance of Air Gap (mm)         
        
        
        # Number of Projections
        self.nProj = 15  
        
        # Angle settings (Degrees)
        self.tubeAngle = 15         # Tube Angle
        
        self.detAngle = 4.2         # Detector Angle

        self.x_offset = 0       # Volume X offset (mm)
        self.y_offset = 0       # Volume Y offset (mm) 
    
    def GE(self):
         
        # Breast voxels density 
        self.nx = 1058    # number of voxels (columns)
        self.ny = 1978    # number of voxels (rows)
        self.nz = 107     # number of voxels (slices)
        
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

        self.x_offset = 0       # Volume X offset (mm)
        self.y_offset = 0       # Volume Y offset (mm) 

    def Siemens(self):

        # Breast voxels density 
        self.nx = 1372    # number of voxels (columns)
        self.ny = 3264    # number of voxels (rows)
        self.nz = 49      # number of voxels (slices)
        
        # Detector panel pixel density
        self.nu = 1376    # number of pixels (columns)
        self.nv = 3264    # number of pixels (rows)
        
        # Single voxel real size (mm)
        self.dx = 0.085 
        self.dy = 0.085
        self.dz = 1
        
        # Single detector real size (mm)
        self.du = 0.085   
        self.dv = 0.085
        
        # X-ray source and detector distances
        self.DSD = 650                          # Distance from source to detector (mm)
        self.DSO = 633                          # Distance from source to the top of object (mm)
        self.DDR = 47                           # Distance from detector to pivot (mm)
        self.DSR = self.DSD - self.DDR          # Distance from source to pivot (mm)
        self.DAG = 17                           # Distance of Air Gap (mm)
        
        # Number of Projections
        self.nProj = 25  
        
        # Angle settings (Degrees)
        self.tubeAngle = 50     # Tube Angle
        
        self.detAngle = 0       # Detector Angle

        self.x_offset = 0       # Volume X offset (mm)
        self.y_offset = 0       # Volume Y offset (mm) 

            
    def SheppLogan(self):
         
        # Breast voxels density 
        self.nx = 150    # number of voxels (columns)
        self.ny = 150    # number of voxels (rows)
        self.nz = 128    # number of voxels (slices)
        
        # Detector panel pixel density
        self.nu = 256    # number of pixels (columns)
        self.nv = 448    # number of pixels (rows)
        
        # Single voxel real size (mm)
        self.dx = 0.112 
        self.dy = 0.112
        self.dz = 1
        
        # Single detector real size (mm)
        self.du = 0.14
        self.dv = 0.14
        
        # X-ray source and detector distances
        self.DSD = 700                           # Distance from source to detector (mm)
        self.DSO = 597                           # Distance from source to the top of object (mm)
        self.DDR = 0                             # Distance from detector to pivot (mm)
        self.DSR = self.DSD - self.DDR           # Distance from source to pivot (mm)
        self.DAG = 25                            # Distance of Air Gap (mm)         
        
        
        # Number of Projections
        self.nProj = 15  
        
        # Angle settings (Degrees)
        self.tubeAngle = 15     # Tube Angle
        
        self.detAngle = 4.2     # Detector Angle

        self.x_offset = 0       # Volume X offset (mm)
        self.y_offset = 0       # Volume Y offset (mm) 

        