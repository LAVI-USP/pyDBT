/*
%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This is the header function
%     ---------------------------------------------------------------------
%     Copyright (C) <2019>  <Rodrigo de Barros Vimieiro>
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
%% 3-D Back-projection Branchless Distance Driven Code (CPU-Multithread)
*/

//typedef double user_dataType;
//typedef float user_dataType;
//#define user_mexdataType mxDOUBLE_CLASS
//#define user_mexdataType mxSINGLE_CLASS

#include <stdio.h>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

// Global variable
int nThreads;

void backprojectionDDb(double* const pVolume,
	double* const pProj,
	double* const pTubeAngle,
	double* const pDetAngle,
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

void linspace(double start, 
			  double end, 
			  int num,
			  double* pLinspaced);

void mapBoudaries(double* pBound,
	const int nElem,
	const double valueLeftBound,
	const double sizeElem,
	const double offset);

void mapDet2Slice(double* const pXmapp,
	double* const pYmapp,
	double tubeX,
	double tubeY,
	double tubeZ,
	double * const pXcoord,
	double * const pYcoord,
	double * const pZcoord,
	double ZSlicecoord,
	const int nXelem,
	const int nYelem);

void bilinear_interpolation(double* pSliceI,
	double* pProj,
	double* pObjX,
	double* pObjY,
	double* pDetmX,
	double* pDetmY,
	const int nPixXMap,
	const int nPixYMap,
	const int nDetXMap,
	const int nDetYMap,
	const int nDetX,
	const int nDetY,
	const unsigned int np);

void differentiation(double* pVolume,
	double* pSliceI,
	double tubeX,
	double rtubeY,
	double rtubeZ,
	double* const pObjX,
	double* const pObjY,
	double* const pObjZ,
	const int nPixX,
	const int nPixY,
	const int nPixXMap,
	const int nPixYMap,
	const double du,
	const double dv,
	const double dx,
	const double dy,
	const double dz,
	const unsigned int nz);
