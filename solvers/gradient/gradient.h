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

#include "mesh.h"
#include "mesh2D.h"
#include "mesh3D.h"

// block size for reduction (hard coded)
#define blockSize 256

typedef struct{

  int dim;
  int elementType; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)

  int Nfields;

  hlong totalElements;
  dlong Nblock;

  dfloat *q, *gradientq;

  dfloat *plotInterp;
  int *plotEToV;

  int frame;

  mesh_t *mesh;

  occa::kernel gradientKernel;
  occa::kernel isoSurfaceKernel;
  occa::memory o_q;
  occa::memory o_gradientq;

  occa::memory o_plotEToV;
  occa::memory o_plotInterp;

  //halo data
  dlong haloBytes;
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  occa::memory o_sendBuffer;
  occa::memory o_recvBuffer;
  occa::memory h_sendBuffer;
  occa::memory h_recvBuffer;
  occa::memory o_haloBuffer;

  dlong haloStressesBytes;
  dfloat *sendStressesBuffer;
  dfloat *recvStressesBuffer;
  occa::memory o_haloStressesBuffer;

}gradient_t;

gradient_t *gradientSetup(mesh_t *mesh, setupAide &options);

void gradientError(mesh_t *mesh, dfloat time);

void gradientReport(gradient_t *gradient, dfloat time, setupAide &options);

void gradientPlotVTU(gradient_t *gradient, int isoNtris, dfloat *isoq, char *fileName);

