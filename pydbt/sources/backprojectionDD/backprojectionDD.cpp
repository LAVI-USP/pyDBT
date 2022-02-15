/*
%% Author: Rodrigo de Barros Vimieiro
% Date: January, 2021
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 backprojectionDD_lib(pProj, pVolume, pGeo)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function reconstruct the 3D volume from projections, based on
%     the Distance-Driven principle. It works by calculating the overlap
%     in X and Y axis of the volume and the detector boundaries.
%     The geometry is for DBT with half cone-beam. All parameters are set
%     in "ParameterSettings" code.
%
%     INPUT:
%
%     - pProj = Pointer to 2D projections for each angle
%     - pGeo = Pointer to parameters of all pGeometry
%
%     OUTPUT:
%
%     - pVolume = Pointer to reconstructed volume.
%
%     Reference: Three-Dimensional Digital Tomosynthesis - Yulia
%     Levakhina (2014), Cap 3.6 and 3.7.
%
%     Original Paper: De Man, Bruno, and Samit Basu. "Distance-driven
%     projection and backprojection in three dimensions." Physics in
%     Medicine & Biology (2004).
%
%     ---------------------------------------------------------------------
%     Copyright (C) <2018>  <Rodrigo de Barros Vimieiro>
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
%% 3-D Distance Driven Back-projection Code
*/

#include "backprojectionDD.h"

extern "C" void backprojectionDD_lib(double* const pProj, 
								  	 double* const pVolume, 
								  	 float* const pGeo){


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

	const int x_offset = (const int)pGeo[18];
	const int y_offset = (const int)pGeo[19];

	linspace(-tubeAngle/2, tubeAngle/2, nProj, pTubeAngle);
	linspace(-detAngle/2, detAngle/2, nProj, pDetAngle);

	//printf("Nx:%d Ny:%d Nz:%d \nNu:%d Nv:%d \nDx:%.2f Dy:%.2f Dz:%.2f \nDu:%.2f Dv:%.2f \nDSD:%.2f DDR:%.2f \n", nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR);

	const int nDetXMap = nDetX + 1;
	const int nDetYMap = nDetY + 1;
	const int nPixXMap = nPixX + 1;
	const int nPixYMap = nPixY + 1;

	// Memory allocation for projections coordinates
	double* const pDetX = (double*)malloc(nDetXMap * sizeof(double));
	double* const pDetY = (double*)malloc(nDetYMap * sizeof(double));
	double* const pDetZ = (double*)malloc(nDetYMap * sizeof(double));
	double* const pObjX = (double*)malloc(nPixXMap * sizeof(double));
	double* const pObjY = (double*)malloc(nPixYMap * sizeof(double));
	double* const pObjZ = (double*)malloc(nSlices * sizeof(double));

	// Memory allocation for mapped coordinates
	double* const pDetmY = (double*)malloc(nDetYMap * sizeof(double));
	double* const pDetmX = (double*)malloc(nDetYMap * nDetXMap * sizeof(double));
	double* const pPixmX = (double*)malloc(nPixXMap * sizeof(double));
	double* const pPixmY = (double*)malloc(nPixYMap * sizeof(double));


	// Memory allocation for rotated detecto coords
	double* const pRdetY = (double*)malloc(nDetYMap * sizeof(double));
	double* const pRdetZ = (double*)malloc(nDetYMap * sizeof(double));

	// Map detector and object boudaries
	mapBoudaries(pDetX, nDetXMap, (double)nDetX, -du, 0.0);

	mapBoudaries(pDetY, nDetYMap, nDetY / 2.0, dv, 0.0);

	mapBoudaries(pDetZ, nDetYMap, 0.0, 0.0, 0.0);

	mapBoudaries(pObjX, nPixXMap, (double)nPixX, -dx, x_offset);

	mapBoudaries(pObjY, nPixYMap, nPixY / 2.0, dy, y_offset);

	mapBoudaries(pObjZ, nSlices, 0.0, dz, DAG + (dz / 2.0));


	// X - ray tube initial position
	double tubeX = 0;
	double tubeY = 0;
	double tubeZ = DSD;

	// Iso - center position
	double isoY = 0;
	double isoZ = DDR;


	// Allocate memory for temp projection variable
	double* const pVolumet = (double*)malloc(nPixY *nPixX * nSlices * sizeof(double));


	// Initiate temp volume with zeros
	for (int z = 0; z < nSlices; z++)
		for (int x = 0; x < nPixX; x++)
			for (int y = 0; y < nPixY; y++)
				pVolumet[(z*nPixX*nPixY) + (x*nPixY) + y] = 0;


	// Allocate memory for flipped projection
	double* const pProjf = (double*)malloc(nDetY *nDetX * nProj * sizeof(double));

	// Flip projection X (Img coord is reverse to Global)
	for (int p = 0; p < nProj; p++)
		for (int x = 0, x_inv = nDetX - 1; x < nDetX; x++, x_inv--)
			for (int y = 0; y < nDetY; y++)
				pProjf[(p*nDetX*nDetY) + (x*nDetY) + y] = pProj[(p*nDetX*nDetY) + (x_inv*nDetY) + y];


	// For each projection
	for (int p = 0; p < nProj; p++) {

		// Get specif tube angle for the projection
		double theta = pTubeAngle[p] * M_PI / 180.0;

		// Get specif detector angle for the projection
		double phi = pDetAngle[p] * M_PI / 180.0;


		// Tubre rotation
		double rtubeY = ((tubeY - isoY)*(double)cos(theta) - (tubeZ - isoZ)*(double)sin(theta)) + isoY;
		double rtubeZ = ((tubeY - isoY)*(double)sin(theta) + (tubeZ - isoZ)*(double)cos(theta)) + isoZ;

		// Detector rotation
		for (int v = 0; v < nDetYMap; v++) {
			pRdetY[v] = ((pDetY[v] - isoY)*(double)cos(phi) - (pDetZ[v] - isoZ)*(double)sin(phi)) + isoY;
			pRdetZ[v] = ((pDetY[v] - isoY)*(double)sin(phi) + (pDetZ[v] - isoZ)*(double)cos(phi)) + isoZ;
		}


		// Map detector onto XY plane(Inside proj loop in case detector rotates)
		mapp2xy(pDetmX, pDetmY, tubeX, rtubeY, rtubeZ, pDetX, pRdetY, pRdetZ, nDetXMap, nDetYMap);


		// Pixel start index and increment
		int detIstart = 0;
		int detIinc = 1;

		// Mapped detector length
		double deltaDetmY = pDetmY[detIstart + detIinc] - pDetmY[detIstart];


		// For each slice
		for (int z = 0; z < nSlices; z++) {

			// Tmp Z coords value for Y direction
			double* const pObjZt = (double*)malloc(nPixYMap * sizeof(double));

			// Tmp Pixel X mapped coords
			double* const pPixmXt = (double*)malloc(nPixYMap * nPixXMap * sizeof(double));

			// Get specif Z coord value for each slice
			for (int k = 0; k < nPixYMap; k++)	pObjZt[k] = pObjZ[z];

			// Map slice onto XY plane
			mapp2xy(pPixmXt, pPixmY, tubeX, rtubeY, rtubeZ, pObjX, pObjY, pObjZt, nPixXMap, nPixYMap);

			// Flip X(Img coord is reverse to Global)
			for (int x = 0, x_inv = nPixXMap - 1; x < nPixXMap; x++, x_inv--)
				pPixmX[x] = pPixmXt[x_inv*nPixYMap];

			// Free temp variables
			free(pObjZt);
			free(pPixmXt);

			// Pixel start index and increment
			int pixIstart = 0;
			int pixIinc = 1;

			// Mapped pixel length
			double deltaPixmX = pPixmX[pixIstart + pixIinc] - pPixmX[pixIstart];
			double deltaPixmY = pPixmY[pixIstart + pixIinc] - pPixmY[pixIstart];

			// Start pixel and detector indices
			int detIndY = detIstart;
			int pixIndY = pixIstart;

			// Case 1
			// Find first detector overlap maped with pixel maped on Y
			if (pDetmY[detIndY] - pPixmY[pixIstart] < -deltaDetmY)
				while (pDetmY[detIndY] - pPixmY[pixIstart] < -deltaDetmY)
					detIndY = detIndY + detIinc;

			else
				// Case 2
				// Find first pixel overlap maped with detector maped on Y
				if (pDetmY[detIstart] - pPixmY[pixIndY] > deltaPixmY)
					while (pDetmY[detIstart] - pPixmY[pixIndY] > deltaPixmY)
						pixIndY = pixIndY + pixIinc;

			double moving_left_boundaryY;

			// Get the left coordinate of the first overlap on Y axis
			if (pDetmY[detIndY] < pPixmY[pixIndY])
				moving_left_boundaryY = pPixmY[pixIndY];
			else
				moving_left_boundaryY = pDetmY[detIndY];


			// Allocate memory for specif row of X map detector coords
			double* const pDetmXrow = (double*)malloc(nDetXMap * sizeof(double));

			double overLapY;

			// Loop over Y intersections
			while ((detIndY < nDetY) && (pixIndY < nPixY)) {

				double alpha = (double)atan((pDetmY[detIndY] + (deltaDetmY / 2) - rtubeY) / rtubeZ);

				// Case A, when you jump to the next detector boundarie but stay
				// in the same pixel
				if (pDetmY[detIndY + 1] <= pPixmY[pixIndY + 1])
					overLapY = (pDetmY[detIndY + 1] - moving_left_boundaryY) / deltaDetmY; // Normalized overlap Calculation

				else
					// Case B, when you jump to the next pixel boundarie but stay
					// in the same detector
					overLapY = (pPixmY[pixIndY + 1] - moving_left_boundaryY) / deltaDetmY; // Normalized overlap Calculation

																						   //										***** X overlap *****
				int detIndX = detIstart;
				int pixIndX = pixIstart;


				// Get row / coll of X flipped, which correspond to that Y overlap det
				for (int x = 0, x_inv = nDetXMap - 1; x < nDetXMap; x++, x_inv--)
					pDetmXrow[x] = pDetmX[(x_inv*nDetYMap) + detIndY];

				// Mapped detecor length on X
				double deltaDetmX = pDetmXrow[detIstart + detIinc] - pDetmXrow[detIstart];

				// Case 1
				// Find first detector overlap maped with pixel maped on X
				if (pDetmXrow[detIndX] - pPixmX[pixIstart] < -deltaDetmX)
					while (pDetmXrow[detIndX] - pPixmX[pixIstart] < -deltaDetmX)
						detIndX = detIndX + detIinc;

				else
					// Case 2
					// Find first pixel overlap maped with detector maped on X
					if (pDetmXrow[detIstart] - pPixmX[pixIndX] > deltaPixmX)
						while (pDetmXrow[detIstart] - pPixmX[pixIndY] > deltaPixmX)
							pixIndX = pixIndX + pixIinc;

				double moving_left_boundaryX;

				// Get the left coordinate of the first overlap on X axis
				if (pDetmXrow[detIndX] < pPixmX[pixIndX])
					moving_left_boundaryX = pPixmX[pixIndX];
				else
					moving_left_boundaryX = pDetmXrow[detIndX];


				// Loop over X intersections
				while ((detIndX < nDetX) && (pixIndX < nPixX)) {

					double gamma = (double)atan((pDetmXrow[detIndX] + (deltaDetmX / 2) - tubeX) / rtubeZ);

					// Case A, when you jump to the next detector boundarie but stay
					// in the same pixel
					if (pDetmXrow[detIndX + 1] <= pPixmX[pixIndX + 1]) {

						double overLapX = (pDetmXrow[detIndX + 1] - moving_left_boundaryX) / deltaDetmX; // Normalized overlap Calculation

						pVolumet[(z*nPixX*nPixY) + (pixIndX*nPixY) + pixIndY] += overLapX * overLapY * pProjf[(p*nDetX*nDetY) + (detIndX *nDetY) + detIndY] * dz / ((double)cos(alpha)*(double)cos(gamma));

						detIndX = detIndX + detIinc;
						moving_left_boundaryX = pDetmXrow[detIndX];
					}
					else {
						// Case B, when you jump to the next pixel boundarie but stay
						// in the same detector

						double overLapX = (pPixmX[pixIndX + 1] - moving_left_boundaryX) / deltaDetmX; // Normalized overlap Calculation

						pVolumet[(z*nPixX*nPixY) + (pixIndX*nPixY) + pixIndY] += overLapX * overLapY * pProjf[(p*nDetX*nDetY) + (detIndX *nDetY) + detIndY] * dz / ((double)cos(alpha)*(double)cos(gamma));

						pixIndX = pixIndX + pixIinc;
						moving_left_boundaryX = pPixmX[pixIndX];

					}

				}
				//										***** Back to Y overlap *****

				// Case A, when you jump to the next detector boundarie but stay
				// in the same pixel
				if (pDetmY[detIndY + 1] <= pPixmY[pixIndY + 1]) {
					detIndY = detIndY + detIinc;
					moving_left_boundaryY = pDetmY[detIndY];
				}
				else {
					// Case B, when you jump to the next pixel boundarie but stay
					// in the same detector
					pixIndY = pixIndY + pixIinc;
					moving_left_boundaryY = pPixmY[pixIndY];
				}

			} // Y Overlap loop

			// Free memory
			free(pDetmXrow);

		} // Loop end slices

	} // Loop end Projections


	// Free memory
	free(pProjf);
	free(pDetX);
	free(pDetY);
	free(pDetZ);
	free(pObjX);
	free(pObjY);
	free(pObjZ);
	free(pDetmY);
	free(pDetmX);
	free(pPixmX);
	free(pPixmY);
	free(pRdetY);
	free(pRdetZ);
	free(pTubeAngle);
	free(pDetAngle);

	// Flip volume back X (Img coord is reverse to Global)
	for (int z = 0; z < nSlices; z++)
		for (int x = 0, x_inv = nPixX - 1; x < nPixX; x++, x_inv--)
			for (int y = 0; y < nPixY; y++)
				pVolume[(z*nPixX*nPixY) + (x*nPixY) + y] = pVolumet[(z*nPixX*nPixY) + (x_inv*nPixY) + y] / (double) nProj;

	free(pVolumet);

	return;
}


// Make boundaries of detector and slices
void mapBoudaries(double* pBound, 
	const int nElem, 
	const double valueLeftBound, 
	const double sizeElem, 
	const double offset){

	for (int k = 0; k < nElem; k++)
		pBound[k] = (k - valueLeftBound) * sizeElem + offset;

	return;
}

// Map on XY plane
void mapp2xy(double* const pXmapp,
	double* const pYmapp,
	double tubeX,
	double tubeY,
	double tubeZ,
	double * const pXcoord,
	double * const pYcoord,
	double * const pZcoord,
	const int nXelem,
	const int nYelem){


	for (int x = 0; x < nXelem; x++)
		for (int y = 0; y < nYelem; y++) {

			int ind = (x*nYelem) + y;
			pXmapp[ind] = pXcoord[x] + pZcoord[y] * (pXcoord[x] - tubeX) / (tubeZ - pZcoord[y]);

			if (x == 0)
				pYmapp[y] = pYcoord[y] + pZcoord[y] * (pYcoord[y] - tubeY) / (tubeZ - pZcoord[y]);
		}


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