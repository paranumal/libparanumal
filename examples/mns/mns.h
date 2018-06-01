#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "mesh3D.h"
#include "elliptic.h"

typedef struct {

  int dim, elementType;

  mesh_t *mesh;
  elliptic_t *uSolver;
  elliptic_t *vSolver;
  elliptic_t *wSolver;
  elliptic_t *pSolver;

  setupAide options;
  setupAide vOptions, pOptions; 	

  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu, Re;
  dfloat ubar, vbar, wbar, pbar;
  int NVfields, NTfields;
  dlong fieldOffset;
  dlong Ntotal;

  int Nblock;

  dfloat dt;          // time step
  dfloat dtMIN;         
  dfloat time;
  int tstep;
  dfloat g0, ig0, lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat startTime;   
  dfloat finalTime;   

  int temporalOrder;
  int ExplicitOrder; 
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     
  int   outputStep; 

  int ARKswitch;
  
  int NiterU, NiterV, NiterW, NiterP;
  
  // Reinitialization
  dfloat hmin;
  dfloat *SPhi, *GPhi; 

  //solver tolerances
  dfloat presTOL, velTOL;

  dfloat idt, inu; // hold some inverses
  
  dfloat *U, *P, *Phi;
  dfloat *NU, *LU, *GP;
  dfloat *GU;   
  dfloat *rhsU, *rhsV, *rhsW, *rhsP, *rhsPhi;   
  dfloat *rkU, *rkP, *PI;
  dfloat *rkNU, *rkLU, *rkGP;
  
  dfloat *Vort, *Div;

  dfloat g[3];      // gravitational Acceleration


  //ARK data
  int Nrk;
  dfloat *rkC;
  dfloat *erkA, *irkA, *prkA;
  dfloat *erkB, *irkB, *prkB;
  dfloat *erkE, *irkE, *prkE;
  int embeddedRKFlag;

  //EXTBDF data
  dfloat *extbdfA, *extbdfB, *extbdfC;
  dfloat *extC;

  int *VmapB, *PmapB;
  occa::memory o_VmapB, o_PmapB;

  //halo data
  dfloat *vSendBuffer;
  dfloat *vRecvBuffer;
  dfloat *pSendBuffer;
  dfloat *pRecvBuffer;
  dfloat * velocityHaloGatherTmp;


  int Nsubsteps;  
  dfloat *Ud, *Ue, *resU, *rhsUd, sdt;
  occa::memory o_Ud, o_Ue, o_resU, o_rhsUd;

  dfloat *cU, *cUd;
  occa::memory o_cU, o_cUd;

  occa::kernel scaledAddKernel;
  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleExtKernel;


  occa::memory o_U, o_P, o_Phi;
  occa::memory o_rhsU, o_rhsV, o_rhsW, o_rhsP; 
  occa::memory o_resPhi, o_rhsPhi, o_SPhi, o_GPhi; 

  occa::memory o_NU, o_LU, o_GP;
  occa::memory o_GU;

  occa::memory o_UH, o_VH, o_WH;
  occa::memory o_rkU, o_rkP, o_PI;
  occa::memory o_rkNU, o_rkLU, o_rkGP;

  occa::memory o_Vort, o_Div;

  occa::memory o_vHaloBuffer, o_pHaloBuffer; 
  occa::memory o_velocityHaloGatherTmp;

  //ARK data
  occa::memory o_rkC;
  occa::memory o_erkA, o_irkA, o_prkA;
  occa::memory o_erkB, o_irkB, o_prkB;
  occa::memory o_erkE, o_irkE, o_prkE;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

  occa::kernel velocityHaloExtractKernel;
  occa::kernel velocityHaloScatterKernel;
  occa::kernel pressureHaloExtractKernel;
  occa::kernel pressureHaloScatterKernel;

  occa::kernel setFlowFieldKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;

  occa::kernel levelSetVolumeKernel;
  occa::kernel levelSetSurfaceKernel;
  occa::kernel levelSetUpdateKernel;

  occa::kernel diffusionKernel;
  occa::kernel diffusionIpdgKernel;
  occa::kernel velocityGradientKernel;

  occa::kernel gradientVolumeKernel;
  occa::kernel gradientSurfaceKernel;

  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;
  
  
  occa::kernel pressureRhsKernel;
  occa::kernel pressureRhsIpdgBCKernel;
  occa::kernel pressureRhsBCKernel;
  occa::kernel pressureAddBCKernel;
  occa::kernel pressurePenaltyKernel;
  occa::kernel pressureUpdateKernel;

  occa::kernel velocityRhsKernel;
  occa::kernel velocityRhsIpdgBCKernel;
  occa::kernel velocityRhsBCKernel;
  occa::kernel velocityAddBCKernel;
  occa::kernel velocityUpdateKernel;  
  
  occa::kernel vorticityKernel;

  occa::kernel regularizedSignumKernel;
  occa::kernel reinitializationVolumeKernel;
  occa::kernel reinitializationSurfaceKernel;
  occa::kernel reinitializationUpdateKernel;

}mns_t;

mns_t *mnsSetup(mesh_t *mesh, setupAide options);


void mnsPlotVTU(mns_t *mns, char *fileNameBase);
void mnsLevelSetRun(mns_t *mns);
void mnsLevelSetStep(mns_t *mns, int tstep, int haloBytes,dfloat * sendBuffer, dfloat *recvBuffer);

void mnsReinitializationRun(mns_t *mns);
void mnsReinitializationStep(mns_t *mns, int tstep, int haloBytes,dfloat * sendBuffer, dfloat *recvBuffer);
void mnsComputeSignumTerm(mns_t *mns, int haloBytes,dfloat * sendBuffer, dfloat *recvBuffer);



void mnsReport(mns_t *mns, dfloat time,  int tstep);
void mnsError(mns_t *mns, dfloat time);

void mnsMakePeriodic(mesh_t *mesh, dfloat xper, dfloat yper);

#if 0
void insRunARK(ins_t *ins);
void insRunEXTBDF(ins_t *ins);


void insAdvection(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_NU);
void insDiffusion(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_LU);
void insGradient (ins_t *ins, dfloat time, occa::memory o_P, occa::memory o_GP);
void insDivergence(ins_t *ins,dfloat time, occa::memory o_U, occa::memory o_DU);
void insSubCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_U, occa::memory o_NU);

void insVelocityRhs  (ins_t *ins, dfloat time, int stage, occa::memory o_rhsU, occa::memory o_rhsV, occa::memory o_rhsW);
void insVelocitySolve(ins_t *ins, dfloat time, int stage, occa::memory o_rhsU, occa::memory o_rhsV, occa::memory o_rhsW, occa::memory o_rkU);
void insVelocityUpdate(ins_t *ins, dfloat time, int stage, occa::memory o_rkGP, occa::memory o_rkU);

void insPressureRhs  (ins_t *ins, dfloat time, int stage);
void insPressureSolve(ins_t *ins, dfloat time, int stage);
void insPressureUpdate(ins_t *ins, dfloat time, int stage, occa::memory o_rkP);
#endif
