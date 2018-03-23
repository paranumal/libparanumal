#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "ellipticQuad2D.h"

#define UXID 0
#define UYID 1
#define PRID 2

typedef struct {

  mesh_t *mesh;
  solver_t *uSolver;
  solver_t *vSolver;
  solver_t *pSolver;

  char *pSolverOptions, *vSolverOptions; 	
  precon_t *precon;

  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu;
  int NVfields, NTfields, Nfields;
  int NtotalDofs, NDofs, NtotalElements; // Total DOFs for Velocity i.e. Nelements + Nelements_halo
  int ExplicitOrder; 

  dfloat dt;          // time step
  dfloat lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat finalTime;   // final time to run acoustics to
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     // Number of history states to store
  int   index;       // Index of current state
  int   errorStep; 
  
  int NiterU, NiterV, NiterP;

  //solver tolerances
  dfloat presTOL, velTOL;


  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 
  dfloat idt, ig0, inu; // hold some inverses
  
  dfloat *U, *V, *P, *NU, *NV;   
  dfloat *rhsU, *rhsV, *rhsP;   
  dfloat *PI,*Px,*Py;
  
  dfloat *Vort, *Div;   

  dfloat g[2];      // gravitational Acceleration

  int *VmapB, *PmapB;
  occa::memory o_VmapB, o_PmapB;

  //halo data
  dfloat *tSendBuffer;
  dfloat *tRecvBuffer;
  dfloat *vSendBuffer;
  dfloat *vRecvBuffer;
  dfloat *pSendBuffer;
  dfloat *pRecvBuffer;

  int Nsubsteps;  
  dfloat *Ud, *Vd, *Ue, *Ve, *resU, *resV, sdt;
  occa::memory o_Ud, o_Vd, o_Ue, o_Ve, o_resU, o_resV;

  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleExtKernel;

 
  occa::memory o_U, o_V, o_P, o_NU, o_NV;
  occa::memory o_rhsU, o_rhsV, o_rhsP; 

  occa::memory o_Px, o_Py;

  occa::memory o_UH, o_VH;
  occa::memory o_PI, o_PIx, o_PIy;

  occa::memory o_Ut, o_Vt, o_Pt, o_WN; 

  occa::memory o_Vort, o_Div; 

  // multiple RHS pressure projection variables
  int maxPresHistory, NpresHistory;
  int Nblock;

  dfloat *presAlpha, *presLocalAlpha; 
  occa::memory o_presAlpha;

  dfloat *presHistory;
  occa::memory o_PIbar, o_APIbar, o_presHistory;

  dfloat *blockReduction;
  occa::memory o_blockReduction;

  occa::kernel multiWeightedInnerProductKernel;
  occa::kernel multiInnerProductKernel;
  occa::kernel multiScaledAddKernel;



  occa::memory o_vHaloBuffer, o_pHaloBuffer, o_tHaloBuffer; 

  occa::kernel scaledAddKernel;

  occa::kernel totalHaloExtractKernel;
  occa::kernel totalHaloScatterKernel;

  occa::kernel velocityHaloExtractKernel;
  occa::kernel velocityHaloScatterKernel;
  occa::kernel pressureHaloExtractKernel;
  occa::kernel pressureHaloScatterKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionUpdateKernel; // SS
  
  occa::kernel poissonUpdateKernel; //SS
  occa::kernel poissonRhsCurlKernel; // SS
  occa::kernel poissonRhsNeumannKernel; // SS
  occa::kernel poissonRhsForcingKernel;
  occa::kernel poissonRhsIpdgBCKernel;
  occa::kernel poissonRhsBCKernel;
  occa::kernel poissonAddBCKernel;
  occa::kernel poissonPenaltyKernel;

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
  occa::kernel helmholtzRhsBCKernel;
  occa::kernel helmholtzAddBCKernel;
  
  occa::kernel updateUpdateKernel;
  occa::kernel vorticityKernel;

}ins_t;

ins_t *insSetupQuad2D(mesh2D *mesh, int i, char *options, 
                  char *velSolverOptions, char *velParAlmondOptions,
                  char *prSolverOptions,  char *prParAlmondOptions,
                  char *bdryHeaderFileName);

void insRunQuad2D(ins_t *solver, char *options);
void insPlotVTUQuad2D(ins_t *solver, char *fileNameBase);
void insReportQuad2D(ins_t *solver, int tstep, char *options);
void insErrorQuad2D(ins_t *solver, dfloat time, char *options);

void insAdvectionStepQuad2D(ins_t *solver, int tstep, char * options);
void insAdvectionSubCycleStepQuad2D(ins_t *solver, int tstep, char * options);
void insHelmholtzStepQuad2D(ins_t *solver, int tstep, char * options);
void insPoissonStepQuad2D(ins_t *solver, int tstep, char * options);
void insUpdateStepQuad2D(ins_t *solver, int tstep, char * options);
