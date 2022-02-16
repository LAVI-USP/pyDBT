/*
%% Author: Rodrigo de Barros Vimieiro
% Date: Feb, 2022
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This is the header function
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
%% 3-D Distance Driven Back-projection header
*/

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

// Includes CUDA
#include <cuda_runtime.h>

// Utilities and timing functions
#include <helper_functions.h>    // includes cuda.h and cuda_runtime_api.h

// CUDA helper functions
#include <helper_cuda.h>         // helper functions for CUDA error check

#define integrateXcoord 1
#define integrateYcoord 0


/*
Reference: TIGRE - https://github.com/CERN/TIGRE
*/
#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                printf("%s \n",msg);\
        } \
} while (0)

void linspace(double start, 
                double end, 
                int num,
                double* pLinspaced);

void backprojectionDDb(double* const h_pVolume,
	double* const h_pProj,
	double* const h_pTubeAngle,
	double* const h_pDetAngle,
	const unsigned int nProj,
	const unsigned int nPixX,
	const unsigned int nPixY,
	const unsigned int nSlices,
	const unsigned int nDetX,
	const unsigned int nDetY,
	const signed int idXProj,
	const double x_offset,
	const double y_offset,
	const double dx,
	const double dy,
	const double dz,
	const double du,
	const double dv,
	const double DSD,
	const double DDR,
	const double DAG);
