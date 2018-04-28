#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "mesh2D.h"

// block size for reduction (hard coded)
#define blockSize 256

typedef struct{

  int Nstresses;
  int Nfields;

  hlong totalElements;
  dlong Nblock;
  
  dfloat *q, *rhsq, *resq;
  
  dfloat *viscousStresses;
  dfloat *Vort;
  dfloat *LIFTT;

  dfloat *rkq, *rkrhsq, *rkerr;
  dfloat *errtmp;
  int frame;
  
  dfloat mu;
  dfloat RT;
  dfloat rbar;
  dfloat ubar;
  dfloat vbar;
      
  
  mesh_t *mesh;

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel cubatureVolumeKernel;
  occa::kernel cubatureSurfaceKernel;
  occa::kernel updateKernel;
  occa::kernel rkStageKernel;
  occa::kernel rkUpdateKernel;
  occa::kernel rkErrorEstimateKernel;

  occa::kernel stressesVolumeKernel;
  occa::kernel stressesSurfaceKernel;

  occa::kernel vorticityKernel;
  
  occa::memory o_LIFTT;
  
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
  occa::memory o_haloBuffer;

  dlong haloStressesBytes;
  dfloat *sendStressesBuffer;
  dfloat *recvStressesBuffer;
  occa::memory o_haloStressesBuffer;

  // DOPRI5 RK data
  int advSwitch;
  int Nrk;
  dfloat ATOL, RTOL;
  dfloat factor1, invfactor1;
  dfloat factor2, invfactor2;
  dfloat exp1, facold,  dtMIN, safe, beta;
  dfloat *rkA, *rkC, *rkE;
  occa::memory o_rkA, o_rkC, o_rkE;

  
}cns_t;

void cnsRunQuad2D(cns_t *cns, setupAide &newOptions);

cns_t *cnsSetupQuad2D(mesh2D *mesh, setupAide &newOptions, char* boundaryHeaderFileName);

void cnsError2D(mesh_t *mesh, dfloat time);

void cnsCavitySolution2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);


void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);

void cnsReportQuad2D(cns_t *cns, dfloat time, setupAide &newOptions);

void cnsPlotVTUQuad2D(cns_t *cns, char *fileName);

void cnsDopriStepQuad2D(cns_t *cns, setupAide &newOptions, const dfloat time);

void cnsLserkStepQuad2D(cns_t *cns, setupAide &newOoptions, const dfloat time);

dfloat cnsDopriEstimateQuad2D(cns_t *cns);
