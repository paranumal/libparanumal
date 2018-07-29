/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

#define WADG 1
#define ASYNC 1

void acousticsSetup2D(mesh2D *mesh);

void acousticsRun2Dbbdg(mesh2D *mesh);
void acousticsOccaRun2Dbbdg(mesh2D *mesh);

void acousticsMRABUpdate2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt);  
void acousticsMRABUpdateTrace2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt); 
void acousticsMRABUpdate2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt); 
void acousticsMRABUpdateTrace2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt);  

void acousticsVolume2Dbbdg(mesh2D *mesh, int lev);
void acousticsSurface2Dbbdg(mesh2D *mesh, int lev, dfloat t);

void acousticsPmlSetup2D(mesh2D *mesh);

void acousticsPmlVolume2Dbbdg(mesh2D *mesh, int lev);
void acousticsPmlSurface2Dbbdg(mesh2D *mesh, int lev, dfloat t);
void acousticsMRABpmlUpdate2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);
void acousticsMRABpmlUpdateTrace2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);
void acousticsMRABpmlUpdate2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);
void acousticsMRABpmlUpdateTrace2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);

void acousticsSourceSetup2D(mesh2D *mesh, occa::kernelInfo &kernelInfo);


void acousticsError2D(mesh2D *mesh, dfloat time);

void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t, dfloat *u, dfloat *v, dfloat *p);
