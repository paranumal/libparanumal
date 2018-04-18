
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh3D.h"
#include <complex.h>

typedef struct {
  
  mesh_t   *mesh;

  char *options;
  
  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu;

  dfloat dt;          // time step
  dfloat lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat finalTime;   // final time to run acoustics to
  iint   NtimeSteps;  // number of time steps 
  iint   Nstages;     // Number of history states to store
  iint   index;       // Index of current state
  iint   errorStep; 

  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 

#if 0
  dfloat *resQ;
  dfloat *rhsQ;
  dfloat *Q;
#endif
  
  iint Nsubsteps;  
  dfloat *rka;
  dfloat *rkb;
  dfloat *rkc;
  //  occa::memory o_Q, o_resQ, o_rhsQ;
  
  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel rkUpdateKernel;

  occa::memory o_haloBuffer;

  occa::kernel haloExtractKernel;
  occa::kernel haloScatterKernel;

  //DOPRI parameters
  dfloat dtmin; //minumum allowed timestep
  dfloat absTol; //absolute error tolerance
  dfloat relTol; //relative error tolerance
  dfloat safety; //safety factor
  dfloat beta;     //error control parameters
  dfloat factor1;
  dfloat factor2;
  dfloat exp1;
  dfloat invfactor1;
  dfloat invfactor2;
  dfloat oldFactor;
  dfloat outputInterval;
  dfloat nextOutputTime;
  dfloat outputNumber;
  dfloat time;
  iint tstep;
  iint allStep;
  iint Nblock;
  iint blockSize;
  iint Ntotal;
  dfloat *errtmp;
}solver_t;

solver_t *advectionSetupQuad3D(mesh_t *mesh,char* mode);

void advectionRunMRSAABQuad3D(solver_t *solver);
void advectionRunLSERKQuad3D(solver_t *solver);
void advectionRunDOPRIQuad3D(solver_t *solver);

void advectionPlotVTUQuad3D(mesh_t *mesh, char *fileNameBase, iint fld);
void advectionPlotVTUQuad3DV2(mesh_t *mesh, char *fileNameBase, iint tstep);

void advectionOccaSetupQuad3D(mesh_t *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);

void meshMRABSetupQuad3D(mesh3D *mesh, dfloat *EToDT, int maxLevels);

void advectionPlotLevels(mesh_t *mesh, char *fileNameBase, iint tstep,dfloat *q);

void advectionPlotNorms(mesh_t *mesh, char *fileNameBase, iint tstep,dfloat *q);

void advectionErrorNormQuad3D(mesh_t *mesh, dfloat t, char *fileBase, int slice);
