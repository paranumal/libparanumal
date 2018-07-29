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
#include "mesh3D.h"

#define USE_BERN 1
#define WADG 1
#define ASYNC 1


void acousticsSetup3D(mesh3D *mesh);

void acousticsRun3Dbbdg(mesh3D *mesh);
void acousticsOccaRun3Dbbdg(mesh3D *mesh);

void acousticsMRABUpdate3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt);
void acousticsMRABUpdateTrace3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt);
void acousticsMRABUpdate3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt);
void acousticsMRABUpdateTrace3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat t, dfloat dt);

void acousticsVolume3Dbbdg(mesh3D *mesh, int lev);
void acousticsSurface3Dbbdg(mesh3D *mesh, int lev, dfloat time);

void acousticsPmlSetup3D(mesh3D *mesh);

void acousticsPmlVolume3Dbbdg(mesh3D *mesh, int lev);
void acousticsPmlSurface3Dbbdg(mesh3D *mesh, int lev, dfloat t);
void acousticsMRABpmlUpdate3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);
void acousticsMRABpmlUpdateTrace3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);
void acousticsMRABpmlUpdate3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);
void acousticsMRABpmlUpdateTrace3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, int lev, dfloat dt);

void acousticsSourceSetup3D(mesh3D *mesh, occa::kernelInfo &kernelInfo);

void acousticsError3D(mesh3D *mesh, dfloat time);

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void acousticsGaussianPulse3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);
