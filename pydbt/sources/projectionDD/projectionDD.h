/*
%% Author: Rodrigo de Barros Vimieiro
% Date: January, 2021
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This is the header function
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
%% 3-D Distance Driven Projection Code
*/

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

void mapBoudaries(double* pBound, 
                  const int nElem, 
                  const double valueLeftBound, 
                  const double sizeElem, 
                  const double offset);

void mapp2xy(double* const pXmapp,
             double* const pYmapp,
             double tubeX,
             double tubeY,
             double tubeZ,
             double* const pXcoord,
             double* const pYcoord,
             double* const pZcoord,
             const int nXelem,
             const int nYelem);

void linspace(double start, 
             double end, 
             int num,
             double* pLinspaced);