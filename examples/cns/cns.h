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

  occa::memory o_Dmatrices;
  occa::memory o_cubDWmatrices;
  occa::memory o_intLIFTT;
  occa::memory o_intInterpT;

  occa::memory o_intsgeo;

  
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

void cnsRun(cns_t *cns, setupAide &newOptions);

cns_t *cnsSetup(mesh_t *mesh, setupAide &newOptions, char* boundaryHeaderFileName);

void cnsError(mesh_t *mesh, dfloat time);

void cnsCavitySolution(dfloat x, dfloat y, dfloat z, dfloat t,
		       dfloat *u, dfloat *v, dfloat *w, dfloat *p);


void cnsGaussianPulse(dfloat x, dfloat y, dfloat z, dfloat t,
		      dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void cnsReport(cns_t *cns, dfloat time, setupAide &newOptions);

void cnsPlotVTU(cns_t *cns, char *fileName);

void cnsDopriStep(cns_t *cns, setupAide &newOptions, const dfloat time);

void cnsLserkStep(cns_t *cns, setupAide &newOoptions, const dfloat time);

dfloat cnsDopriEstimate(cns_t *cns);

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12
