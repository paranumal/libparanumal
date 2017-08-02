
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh3D.h"
#include "ellipticTet3D.h"

typedef struct {

  mesh_t *mesh;
  solver_t *vSolver;
  solver_t *pSolver;

  char *pSolverOptions, *vSolverOptions; 	
  precon_t *precon;

  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu;
  iint NVfields, NTfields, Nfields;
  iint NtotalDofs, NDofs; // Total DOFs for Velocity i.e. Nelements + Nelements_halo
  iint ExplicitOrder; 

  dfloat dt;          // time step
  dfloat lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat finalTime;   // final time to run acoustics to
  iint   NtimeSteps;  // number of time steps 
  iint   Nstages;     // Number of history states to store
  iint   index;       // Index of current state
  iint   errorStep; 


  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 
  dfloat *rhsU, *rhsV, *rhsW, *rhsP;
  dfloat *U, *V, *W, *P; 

  dfloat dtfactor;

  dfloat *NU, *NV, *NW;
  dfloat *Px, *Py, *Pz;
  dfloat *PI;

  dfloat g[3];      // gravitational Acceleration
 
  // iint Nsubsteps;  
  // dfloat *Ud, *Vd, *Ue, *Ve, *resU, *resV, sdt;
  // occa::memory o_Ud, o_Vd, o_Ue, o_Ve, o_resU, o_resV;

  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleExtKernel;

 
  occa::memory o_U, o_V, o_W, o_P;
  occa::memory o_rhsU, o_rhsV, o_rhsW, o_rhsP; 

  occa::memory o_NU, o_NV, o_NW;
  occa::memory o_Px, o_Py, o_Pz;

  occa::memory o_UH, o_VH, o_WH;
  occa::memory o_PI, o_PIx, o_PIy, o_PIz;

  occa::memory o_vHaloBuffer, o_pHaloBuffer, o_tHaloBuffer; 

  occa::kernel totalHaloExtractKernel;
  occa::kernel totalHaloScatterKernel;

  occa::kernel velocityHaloExtractKernel;
  occa::kernel velocityHaloScatterKernel;
  occa::kernel pressureHaloExtractKernel;
  occa::kernel pressureHaloScatterKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;

  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;
  //
  occa::kernel gradientVolumeKernel;
  occa::kernel gradientSurfaceKernel;
  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;
  //
  occa::kernel helmholtzRhsForcingKernel;
  occa::kernel helmholtzRhsIpdgBCKernel;

  occa::kernel poissonRhsForcingKernel;
  occa::kernel poissonRhsIpdgBCKernel;
  occa::kernel poissonPenaltyKernel;
  
  occa::kernel updateUpdateKernel;

}ins_t;

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

} nonZero_t;



ins_t *insSetup3D(mesh3D *mesh,char *options, char *velSolverOptions, char *prSolverOptions, char *bdryHeaderFileName);

void insRun3D(ins_t *solver, char *options);
void insReport3D(ins_t *solver, iint tstep, char *options);
void insError3D(ins_t *solver, dfloat time, char *options);
void insPlotVTU3D(ins_t *solver, char *fileNameBase);
// void insErrorNorms2D(ins_t *solver, dfloat time, char *options);

void insAdvectionStep3D(ins_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);

void insHelmholtzStep3D(ins_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);

void insPoissonStep3D(ins_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);


void insUpdateStep3D(ins_t *solver, iint tstep, iint haloBytes,
	                   dfloat * sendBuffer, dfloat *recvBuffer, char * options);


void insErrorNorms3D(ins_t *solver, dfloat time, char *options);

// void insAdvectionSubCycleStep2D(ins_t *solver, iint tstep,
//                      dfloat * tsendBuffer, dfloat *trecvBuffer, 
//                      dfloat * sendBuffer, dfloat *recvBuffer,char * options);


