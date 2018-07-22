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

  occa::memory o_q;
  occa::memory o_rhsq;
  occa::memory o_resq;
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

  // DOPRI5 RK data
  int advSwitch;
  int Nrk;
  dfloat ATOL, RTOL;
  dfloat factor1, invfactor1;
  dfloat factor2, invfactor2;
  dfloat exp1, facold,  dtMIN, safe, beta;
  dfloat *rkA, *rkC, *rkE;
  occa::memory o_rkA, o_rkC, o_rkE;
  
}acoustics_t;

void acousticsRun(acoustics_t *acoustics, setupAide &newOptions);

acoustics_t *acousticsSetup(mesh_t *mesh, setupAide &newOptions, char* boundaryHeaderFileName);

void acousticsError(mesh_t *mesh, dfloat time);

void acousticsCavitySolution(dfloat x, dfloat y, dfloat z, dfloat t,
		       dfloat *u, dfloat *v, dfloat *w, dfloat *p);


void acousticsGaussianPulse(dfloat x, dfloat y, dfloat z, dfloat t,
		      dfloat *u, dfloat *v, dfloat *w, dfloat *p);

void acousticsReport(acoustics_t *acoustics, dfloat time, setupAide &newOptions);

void acousticsPlotVTU(acoustics_t *acoustics, char *fileName);

void acousticsDopriStep(acoustics_t *acoustics, setupAide &newOptions, const dfloat time);

void acousticsLserkStep(acoustics_t *acoustics, setupAide &newOoptions, const dfloat time);

dfloat acousticsDopriEstimate(acoustics_t *acoustics);

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12
