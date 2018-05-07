#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "mesh3D.h"
#include "elliptic.h"

typedef struct {

  mesh_t *mesh;
  elliptic_t *uSolver;
  elliptic_t *vSolver;
  elliptic_t *wSolver;
  elliptic_t *pSolver;

  setupAide pOptions, vOptions; 	

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
  
  int NiterU, NiterV, NiterW, NiterP;

  //solver tolerances
  dfloat presTOL, velTOL;

  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 
  dfloat idt, ig0, inu; // hold some inverses
  
  dfloat *U, *V, *W, *P, *NU, *NV, *NW;   
  dfloat *rhsU, *rhsV, *rhsW, *rhsP;   
  dfloat *PI,*Px,*Py,*Pz;
  
  dfloat *Vort, *Div;
  dfloat *Vx,*Vy,*Vz;

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

  int Nsubsteps;  
  dfloat *Ud, *Vd, *Wd, *Ue, *Ve, *We, *resU, *resV, *resW, sdt;
  occa::memory o_Ud, o_Vd, o_Wd, o_Ue, o_Ve, o_We, o_resU, o_resV, o_resW;

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

  occa::memory o_Vort, o_Div;
  occa::memory o_Vx, o_Vy, o_Vz;

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

ins_t *insSetup(mesh_t *mesh, int i, setupAide options);

void insRun(ins_t *ins);
void insPlotVTU(ins_t *ins, char *fileNameBase);
void insReport(ins_t *ins, int tstep);
void insError(ins_t *ins, dfloat time);

void insAdvectionStep(ins_t *ins, int tstep);
void insAdvectionSubCycleStep(ins_t *ins, int tstep);
void insHelmholtzStep(ins_t *ins, int tstep);
void insPoissonStep(ins_t *ins, int tstep);
void insUpdateStep(ins_t *ins, int tstep);
