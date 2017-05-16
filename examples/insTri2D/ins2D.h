
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
iint Nfields;
iint NtotalDofs; // Total DOFs for Velocity i.e. Nelements + Nelements_halo
//
dfloat dt; // time step
dfloat finalTime; // final time to run acoustics to
iint   NtimeSteps;// number of time steps 
iint   errorStep; 

dfloat a0, a1, b0, b1, g0; 
dfloat *U,  *Pr, *rhsU, *rhsPr;
dfloat *UO, *UI, *NU; 
dfloat g[2]; // gravitational Acceleration

occa::memory o_U, o_UO, o_UI, o_NU, o_Pr;
occa::memory o_rhsU, o_rhsPr; 


occa::kernel haloExtractKernel;
occa::kernel haloScatterKernel;

//
occa::kernel advectionVolumeKernel; // deprecated
occa::kernel advectionSurfaceKernel; 
occa::kernel advectionUpdateKernel;  
//
occa::kernel pressureRhsVolumeKernel; 
occa::kernel pressureRhsSurfaceKernel;

  
}solver_t;



solver_t *insSetup2D(mesh2D *mesh, char *options);

void insMakePeriodic2D(mesh2D *mesh, dfloat xper, dfloat yper);

void insRun2D(solver_t *solver, char *options);
void insPlotVTU2D(solver_t *solver, char *fileNameBase);
void insReport2D(solver_t *solver, iint tstep, char *options);
void insError2D(solver_t *solver, dfloat time,char *options);

void insAdvectionStep2D(solver_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);

void insPressureStep2D(solver_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);










