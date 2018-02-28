
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh3D.h"

typedef struct {
  
  mesh_t   *mesh;

  char *options;
  
  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu;

  dfloat dt;          // time step
  dfloat lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat finalTime;   // final time to run acoustics to
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     // Number of history states to store
  int   index;       // Index of current state
  int   errorStep; 

  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 

#if 0
  dfloat *resQ;
  dfloat *rhsQ;
  dfloat *Q;
#endif
  
  int Nsubsteps;  

  //  occa::memory o_Q, o_resQ, o_rhsQ;
  
  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel rkUpdateKernel;

  occa::memory o_haloBuffer;

  occa::kernel haloExtractKernel;
  occa::kernel haloScatterKernel;

}solver_t;


solver_t *boltzmannSetupTri3D(mesh_t *mesh);

void boltzmannRunTri3D(solver_t *solver);

void boltzmannPlotVTUTri3D(mesh_t *mesh, char *fileNameBase, int fld);
void boltzmannPlotVTUTri3DV2(mesh_t *mesh, char *fileNameBase, int tstep);

void boltzmannOccaSetupTri3D(mesh_t *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);

