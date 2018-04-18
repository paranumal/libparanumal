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
  
  occa::memory o_q;
  occa::memory o_rhsq;
  occa::memory o_resq;
  occa::memory o_Vort;
  occa::memory o_viscousStresses;

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
  
}cns_t;

void cnsRunTri2D(cns_t *cns, char *options);

cns_t *cnsSetupTri2D(mesh2D *mesh, char* options, char* boundaryHeaderFileName);

void cnsError2D(mesh2D *mesh, dfloat time);

void cnsCavitySolution2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);


void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);

void cnsReportTri2D(cns_t *cns, dfloat time, char *options);

void cnsPlotVTUTri2D(cns_t *cns, char *fileName);
