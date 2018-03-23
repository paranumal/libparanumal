
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh3D.h"
#include "ellipticTet3D.h"

typedef struct {

  mesh_t *mesh;
  solver_t *uSolver;
  solver_t *vSolver;
  solver_t *wSolver;
  solver_t *pSolver;

  char *pSolverOptions, *vSolverOptions; 	
  precon_t *precon;

  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu;
  int NVfields, NTfields, Nfields;
  int NtotalDofs, NDofs; // Total DOFs for Velocity i.e. Nelements + Nelements_halo
  int ExplicitOrder; 

  dfloat dt;          // time step
  dfloat lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat finalTime;   // final time to run acoustics to
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     // Number of history states to store
  int   index;       // Index of current state
  int   errorStep; 
  int   Nsubsteps;  
  //solver tolerances
  int NiterU, NiterV, NiterW, NiterP;

  dfloat presTOL, velTOL;

  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 
  dfloat idt, ig0, inu; // hold some inverses

  dfloat *rhsU, *rhsV, *rhsW, *rhsP;
  dfloat *U, *V, *W, *P; 
  
  dfloat *Ud, *Vd, *Wd, *Ue, *Ve, *We, *resU, *resV, *resW, sdt;
  occa::memory o_Ud, o_Vd, o_Wd, o_Ue, o_Ve, o_We,  o_resU, o_resV, o_resW;

  dfloat *NU, *NV, *NW;
  dfloat *Px, *Py, *Pz;
  dfloat *PI;

  dfloat *Vx, *Vy, *Vz, *Div; 

  dfloat g[3];      // gravitational Acceleration
 
  int *VmapB, *PmapB;
  occa::memory o_VmapB, o_PmapB;

  //halo data
  dfloat *tSendBuffer;
  dfloat *tRecvBuffer;
  dfloat *vSendBuffer;
  dfloat *vRecvBuffer;
  dfloat *pSendBuffer;
  dfloat *pRecvBuffer;

  // int Nsubsteps;  
  // dfloat *Ud, *Vd, *Ue, *Ve, *resU, *resV, sdt;
  // occa::memory o_Ud, o_Vd, o_Ue, o_Ve, o_resU, o_resV;

  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleExtKernel;

 
  occa::memory o_U, o_V, o_W, o_P;
  occa::memory o_rhsU, o_rhsV, o_rhsW, o_rhsP; 

  occa::memory o_NU, o_NV, o_NW;
  occa::memory o_Px, o_Py, o_Pz;

  occa::memory o_UH, o_VH, o_WH;
  occa::memory o_PI, o_PIx, o_PIy, o_PIz;

  occa::memory o_Vx, o_Vy, o_Vz, o_Div;

  int Nblock;

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

  occa::kernel poissonRhsForcingKernel;
  occa::kernel poissonRhsIpdgBCKernel;
  occa::kernel poissonRhsBCKernel;
  occa::kernel poissonAddBCKernel;
  occa::kernel poissonPenaltyKernel;
  
  occa::kernel updateUpdateKernel;
  occa::kernel vorticityKernel;

}ins_t;


ins_t *insSetupTet3D(mesh3D *mesh, int Ns, char *options, 
                  char *velSolverOptions, char *velParAlmondOptions, 
                  char *prSolverOptions,  char *prParAlmondOptions,
                  char *bdryHeaderFileName);

void insRunTet3D(ins_t *solver, char *options);
void insReportTet3D(ins_t *solver, int tstep, char *options);
void insErrorTet3D(ins_t *solver, dfloat time, char *options);
void insErrorNormsTet3D(ins_t *solver, dfloat time, char *options);
void insPlotVTUTet3D(ins_t *solver, char *fileNameBase);
void insPlotSliceTet3D(ins_t *ins, char *fileName, const int Nslices, const char** dim, const dfloat* c);
void insPlotContourTet3D(ins_t *ins, char *fileName, const char*options);

void insAdvectionStepTet3D(ins_t *solver, int tstep, const char * options);
void insAdvectionSubCycleStepTet3D(ins_t *solver, int tstep, const char * options);
void insHelmholtzStepTet3D(ins_t *solver, int tstep, const char * options);
void insPoissonStepTet3D(ins_t *solver, int tstep, const char * options);
void insUpdateStepTet3D(ins_t *solver, int tstep, const char * options);
