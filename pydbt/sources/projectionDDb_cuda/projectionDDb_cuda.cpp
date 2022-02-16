/*
%% Author: Rodrigo de Barros Vimieiro
% Date: Feb, 2022
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projectionDDb_cuda_lib(pVolume, pProj, pGeo, idXProj)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function calculates the volume projection based on the
%     Branchless Distance-Driven principle.
%     The geometry is for DBT with half cone-beam. All parameters are set
%     in "ParameterSettings" code.
%
%     INPUT:
%
%     - pVolume = 3D volume for projection
%     - pGeo = Parameter of all geometry
%	  - idXProj = projection number to be projected
%
%     OUTPUT:
%
%     - pProj = projections for each angle.
%
%     Reference:
%     - Branchless Distance Driven Projection and Backprojection,
%     Samit Basu and Bruno De Man (2006)
%     - GPU Acceleration of Branchless Distance Driven Projection and
%     Backprojection, Liu et al (2016)
%     - GPU-Based Branchless Distance-Driven Projection and Backprojection,
%     Liu et al (2017)
%     - A GPU Implementation of Distance-Driven Computed Tomography,
%     Ryan D. Wagner (2017)
%     ---------------------------------------------------------------------
%     Copyright (C) <2022>  <Rodrigo de Barros Vimieiro>
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
% =========================================================================
%% 3-D Projection Branchless Distance Driven Code (CUDA)
*/

#include "projectionDDb_cuda.hpp"


extern "C" void projectionDDb_cuda_lib(double* const pVolume, 
									   double* const pProj, 
									   float* const pGeo,
									   const signed int idXProj){


	const int unsigned nPixX = (const int)pGeo[0]; 
	const int unsigned nPixY = (const int)pGeo[1]; 
	const int unsigned nSlices = (const int)pGeo[2]; 
	const int unsigned nDetX = (const int)pGeo[3]; 
	const int unsigned nDetY = (const int)pGeo[4]; 

	const double dx = (const double)pGeo[5]; 
	const double dy = (const double)pGeo[6]; 
	const double dz = (const double)pGeo[7]; 
	const double du = (const double)pGeo[8]; 
	const double dv = (const double)pGeo[9]; 

	const double DSD = (const double)pGeo[10]; 
	const double DDR = (const double)pGeo[12]; 
	const double DAG = (const double)pGeo[14]; 

	const int unsigned nProj = (const int)pGeo[15];

	const double tubeAngle = (const double)pGeo[16];
	const double detAngle = (const double)pGeo[17];

	double* const pTubeAngle = (double*)malloc(nProj * sizeof(double));
	double* const pDetAngle = (double*)malloc(nProj * sizeof(double));

	const double x_offset = (const double)pGeo[18];
	const double y_offset = (const double)pGeo[19];

	linspace(-tubeAngle/2, tubeAngle/2, nProj, pTubeAngle);
	linspace(-detAngle/2, detAngle/2, nProj, pDetAngle);

	// printf("Nx:%d Ny:%d Nz:%d \nNu:%d Nv:%d \nDx:%.2f Dy:%.2f Dz:%.2f \nDu:%.2f Dv:%.2f \nDSD:%.2f DDR:%.2f \nTube angle:%.2f \nDet angle:%.2f", nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR, tubeAngle, detAngle);

	projectionDDb(pProj, pVolume, pTubeAngle, pDetAngle, nProj, nPixX, nPixY, nSlices, nDetX, nDetY, idXProj, x_offset, y_offset, dx, dy, dz, du, dv, DSD, DDR, DAG);

	free(pTubeAngle);
	free(pDetAngle);

	return;

}

// Linear spaced vector
// Ref: https://stackoverflow.com/a/27030598/8682939
void linspace(double start, 
			  double end, 
			  int num,
			  double* pLinspaced){

	double delta;

	if(num == 0)
		delta = 0;
	else{
		if(num == 1)
			delta = 1;
		else
			delta = (end - start) / (num - 1);
	}

	if((abs(start) < 0.00001) && (abs(end) < 0.00001))
		delta = 0;

	for (int k = 0; k < num; k++){
			pLinspaced[k] = start + k * delta;
	}

  return;
}