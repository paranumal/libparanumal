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

  dfloat *q, *rhsq, *resq;

  dfloat *Vort;

  dfloat *rkq, *rkrhsq, *rkerr;
  dfloat *errtmp;
  int frame;

  mesh_t *mesh;

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel updateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkUpdateKernel;
  occa::kernel rkErrorEstimateKernel;

  occa::kernel combinedKernel;

  occa::kernel invertMassMatrixKernel;
  occa::kernel invertMassMatrixCombinedKernel;
  
  occa::memory o_q;
  occa::memory o_qtmp0;
  occa::memory o_qtmp1;
  occa::memory o_qtmp2;
  occa::memory o_rhsq;
  occa::memory o_resq;
  occa::memory o_saveq;
  
  // [J*W*c_x, J*W*c_y, J*W*c_z]
  occa::memory o_advectionVelocityJW;

  occa::memory o_cubAdvectionVelocityJW;

  // [Jsurf*Wsurf/(Jvol*Wvol)*(c.n + |c.n|)/2
  occa::memory o_advectionVelocityM;

  // [Jsurf*Wsurf/(Jvol*Wvol)*(c.n - |c.n|)/2
  occa::memory o_advectionVelocityP;

  occa::memory o_diagInvMassMatrix;
  
  occa::memory o_rkq, o_rkrhsq, o_rkerr;
  occa::memory o_errtmp;
  
  //halo data
  dlong haloBytes;
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  occa::memory o_sendBuffer;
  occa::memory o_recvBuffer;
  occa::memory o_haloBuffer;

  // DOPRI5 RK data
  int advSwitch;
  int Nrk;
  dfloat ATOL, RTOL;
  dfloat factor1, invfactor1;
  dfloat factor2, invfactor2;
  dfloat exp1, facold,  dtMIN, safe, beta;
  dfloat *rkA, *rkC, *rkE;
  occa::memory o_rkA, o_rkC, o_rkE;
  
}advection_t;

void advectionRun(advection_t *advection, setupAide &newOptions);

advection_t *advectionSetup(mesh_t *mesh, setupAide &newOptions, char* boundaryHeaderFileName);

void advectionError(advection_t *advection, dfloat time);

void advectionCavitySolution(dfloat x, dfloat y, dfloat z, dfloat t, dfloat *q);

void advectionGaussianPulse(dfloat x, dfloat y, dfloat z, dfloat t, dfloat *q);

void advectionReport(advection_t *advection, dfloat time, setupAide &newOptions);

void advectionPlotVTU(advection_t *advection, char *fileName);

void advectionDopriStep(advection_t *advection, setupAide &newOptions, const dfloat time);

void advectionLserkStep(advection_t *advection, setupAide &newOoptions, const dfloat time);

dfloat advectionDopriEstimate(advection_t *advection);

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12
