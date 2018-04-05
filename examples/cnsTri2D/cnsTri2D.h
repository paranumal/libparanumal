#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "mesh2D.h"

typedef struct{

  int Nstresses;
  int Nfields;

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

  occa::kernel stressesVolumeKernel;
  occa::kernel stressesSurfaceKernel;
  
  occa::kernel vorticityKernel;

  dfloat *viscousStresses;
  dfloat *Vort;
  
  occa::memory o_q;
  occa::memory o_rhsq;
  occa::memory o_resq;
  occa::memory o_Vort;
  occa::memory o_viscousStresses;

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

void cnsReportTri2D(cns_t *cns, int tstep, char *options);

void cnsPlotVTUTri2D(cns_t *cns, char *fileName);
