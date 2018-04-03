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
  occa::kernel updateKernel;

  occa::kernel stressesVolumeKernel;
  occa::kernel stressesSurfaceKernel;

  dfloat *viscousStresses;
  dfloat *LIFTT;
  
  occa::memory o_LIFTT;
  
  occa::memory o_q;
  occa::memory o_rhsq;
  occa::memory o_resq;
  occa::memory o_viscousStresses;

  occa::memory o_haloStressesBuffer;
  
}cns_t;

void cnsRunQuad2D(cns_t *cns);

cns_t *cnsSetupQuad2D(mesh2D *mesh);

void cnsError2D(mesh_t *mesh, dfloat time);

void cnsCavitySolution2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);


void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);

void cnsPlotVTUQuad2D(cns_t *cns, char *fileName);
