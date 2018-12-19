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
  
  int Nstresses;
  int Nfields;

  hlong totalElements;
  dlong Nblock;

  dfloat *q, *rhsq, *resq;

  dfloat *viscousStresses;
  dfloat *Vort;

  dfloat *rkq, *rkrhsq, *rkerr;
  dfloat *errtmp;
  int frame;

  dfloat mu;
  dfloat RT;
  dfloat rbar;
  dfloat ubar;
  dfloat vbar;
  dfloat wbar;

  int outputForceStep;
  
  
  mesh_t *mesh;

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel cubatureVolumeKernel;
  occa::kernel cubatureSurfaceKernel;
  occa::kernel updateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkUpdateKernel;
  occa::kernel rkOutputKernel;
  occa::kernel rkErrorEstimateKernel;

  occa::kernel stressesVolumeKernel;
  occa::kernel stressesSurfaceKernel;
  
  occa::kernel vorticityKernel;

  occa::kernel constrainKernel;
  
  occa::memory o_q;
  occa::memory o_rhsq;
  occa::memory o_resq;
  occa::memory o_Vort;
  occa::memory o_viscousStresses;
  occa::memory o_saveq;
  
  occa::memory o_rkq, o_rkrhsq, o_rkerr;
  occa::memory o_errtmp;

  
  //halo data
  dlong haloBytes;
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  occa::memory o_sendBuffer;
  occa::memory o_recvBuffer;
  occa::memory o_haloBuffer;

  dlong haloStressesBytes;
  dfloat *sendStressesBuffer;
  dfloat *recvStressesBuffer;
  occa::memory o_sendStressesBuffer;
  occa::memory o_recvStressesBuffer;
  occa::memory o_haloStressesBuffer;

  // DOPRI5 RK data
  int advSwitch;
  int Nrk;
  dfloat ATOL, RTOL;
  dfloat factor1, invfactor1;
  dfloat factor2, invfactor2;
  dfloat exp1, facold,  dtMIN, dtMAX, safe, beta;
  dfloat *rkA, *rkC, *rkE, *rkoutB;
  occa::memory o_rkA, o_rkC, o_rkE, o_rkoutB;
  
}cns_t;

void cnsRun(cns_t *cns, setupAide &options);

cns_t *cnsSetup(mesh_t *mesh, setupAide &options);

void cnsError(cns_t *cns, dfloat time);
void cnsForces(cns_t *cns, dfloat time);

void cnsCavitySolution(dfloat x, dfloat y, dfloat z, dfloat t,
		       dfloat *u, dfloat *v, dfloat *w, dfloat *p);


void cnsGaussianPulse(dfloat x, dfloat y, dfloat z, dfloat t,
		      dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void cnsReport(cns_t *cns, dfloat time, setupAide &options);

void cnsPlotVTU(cns_t *cns, char *fileName);

void cnsDopriStep(cns_t *cns, setupAide &options, const dfloat time);
void cnsDopriOutputStep(cns_t *cns, const dfloat time, const dfloat dt, const dfloat outTime, occa::memory o_outq);

void cnsLserkStep(cns_t *cns, setupAide &newOoptions, const dfloat time);

dfloat cnsDopriEstimate(cns_t *cns);

void cnsBodyForce(dfloat t, dfloat *fx, dfloat *fy, dfloat *fz,
		  dfloat *intfx, dfloat *intfy, dfloat *intfz);

void cnsBrownMinionQuad3D(cns_t *cns);

#ifdef RENDER

void simpleRayTracer(int     plotNelements,
		     dfloat *plotx,
		     dfloat *ploty,
		     dfloat *plotz,
		     dfloat *plotq,
		     const char	 *fileBaseName,
		      const int fileIndex);

void cnsRenderQuad3D(cns_t *cns, char *fileBaseName, int fileIndex);

#endif

