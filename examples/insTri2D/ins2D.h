
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"

#define UXID 0
#define UYID 1

typedef struct {

mesh_t *mesh;	

// INS SOLVER OCCA VARIABLES
dfloat rho, nu ;
iint NVfields, NTfields;
iint NtotalDofs, NDofs; // Total DOFs for Velocity i.e. Nelements + Nelements_halo
//
dfloat dt; // time step
dfloat finalTime; // final time to run acoustics to
iint   NtimeSteps;// number of time steps 
iint   errorStep; 

dfloat a0, a1, a2, b0, b1, b2, g0, tau; 
dfloat *Pr, *rhsU, *rhsPr;
dfloat *U, *UO, *UOO, *UI, *NU, *NUO,*NUOO; 
dfloat g[2]; // gravitational Acceleration

occa::memory o_U, o_UO, o_UOO, o_UI, o_NU, o_NUO, o_NUOO, o_Pr;
occa::memory o_rhsU, o_rhsPr; 


occa::kernel helmholtzHaloExtractKernel;
occa::kernel helmholtzHaloScatterKernel;
occa::kernel poissonHaloExtractKernel;
occa::kernel poissonHaloScatterKernel;
//
occa::kernel helmholtzRhsVolumeKernel;
occa::kernel helmholtzRhsSurfaceKernel;
occa::kernel helmholtzRhsUpdateKernel;
occa::kernel helmholtzRhsIpdgBCKernel;

occa::kernel poissonRhsVolumeKernel;
occa::kernel poissonRhsSurfaceKernel;
occa::kernel poissonRhsIpdgBCKernel;



// //
// occa::kernel advectionVolumeKernel; // deprecated
// occa::kernel advectionSurfaceKernel; 
// occa::kernel advectionUpdateKernel;  
// //
// occa::kernel pressureRhsVolumeKernel; 
// occa::kernel pressureRhsSurfaceKernel;

  
}solver_t;



solver_t *insSetup2D(mesh2D *mesh, char *options);

void insMakePeriodic2D(mesh2D *mesh, dfloat xper, dfloat yper);

void insRun2D(solver_t *solver, char *options);
void insPlotVTU2D(solver_t *solver, char *fileNameBase);
void insReport2D(solver_t *solver, iint tstep, char *options);
void insError2D(solver_t *solver, dfloat time,char *options);

void insHelmholtzStep2D(solver_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);

void insPoissonStep2D(solver_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);










