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

