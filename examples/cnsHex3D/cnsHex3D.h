#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "mesh3D.h"

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
  dfloat wbar;
  
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

void cnsRunHex3D(cns_t *cns, setupAide &newOptions);

cns_t *cnsSetupHex3D(mesh3D *mesh, setupAide &newOptions, char* boundaryHeaderFileName);

void cnsError3D(mesh_t *mesh, dfloat time);

void cnsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat t,
			 dfloat *u, dfloat *v, dfloat *w, dfloat *p);


void cnsGaussianPulse3D(dfloat x, dfloat y, dfloat z, dfloat t,
			dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void cnsReportHex3D(cns_t *cns, dfloat time, setupAide &newOptions);

void cnsPlotVTUHex3D(cns_t *cns, char *fileName);

void cnsDopriStepHex3D(cns_t *cns, setupAide &newOptions, const dfloat time);

void cnsLserkStepHex3D(cns_t *cns, setupAide &newOoptions, const dfloat time);

dfloat cnsDopriEstimateHex3D(cns_t *cns);
