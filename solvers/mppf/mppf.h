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
  elliptic_t *phiSolver;
  elliptic_t *psiSolver;

  setupAide options;
  setupAide vOptions, pOptions, phiOptions; // psiOptions; 	

  // INS SOLVER OCCA VARIABLES
  dfloat rho1, mu1, rho2, mu2, rho0, nu0, Re;
  dfloat ubar, vbar, wbar, pbar;
  int NVfields, NTfields, NPfields; // # of velocity, velocity + pressure and Phi fields
  dlong fieldOffset;
  dlong Ntotal;

  int Nblock;

  dfloat dt, cfl, dti;          // time step
  dfloat dtMIN;         
  dfloat time;
  int tstep, frame;
  dfloat g0, ig0, lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat startTime;   
  dfloat finalTime; 

  // Cahn-Hilliard
  dfloat chM, chL, chS, chA, chSeta2; // Mobility, energy Density (lamda), derived parameters (S and alpha)
  dfloat lambdaVel, lambdaPhi, lambdaPsi;   

  int temporalOrder;
  int ExplicitOrder; 
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     
  int   outputStep;
  int   outputForceStep; 
  int   dtAdaptStep; 


  int ARKswitch;
  
  int NiterU, NiterV, NiterW, NiterP, NiterPhi, NiterPsi;


  //solver tolerances
  dfloat presTOL, velTOL, phiTOL;

  dfloat idt, inu; // hold some inverses
  dfloat hmin, eta, eta2, inveta2, factorS;     // element length, interface thickness, eta^2. 
  dfloat *U, *P, *Phi, *Psi, *Rho, *Mu, *GMu;
  dfloat *NU, *LU, *GP, *NPhi, *HPhi, *GHPhi, *GPhi;
  dfloat *GU, *DU, *SU;   // GU is depreceated
  dfloat *rhsU, *rhsV, *rhsW, *rhsP, *rhsPhi;   
  dfloat *rkU, *rkP, *rkPhi, *PI;
  dfloat *rkNU, *rkLU, *rkGP;
  
  dfloat *Vort, *Div, *Uhat;

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
  occa::memory o_VmapB, o_PmapB, o_FmapB;

  //halo data
  dfloat *vSendBuffer;
  dfloat *vRecvBuffer;
  dfloat *pSendBuffer;
  dfloat *pRecvBuffer;
  dfloat *phiSendBuffer;
  dfloat *phiRecvBuffer;
  dfloat *gPhiSendBuffer;
  dfloat *gPhiRecvBuffer;
  dfloat * velocityHaloGatherTmp;

  occa::memory o_vSendBuffer;
  occa::memory o_vRecvBuffer;
  occa::memory o_pSendBuffer;
  occa::memory o_pRecvBuffer;
  occa::memory o_phiSendBuffer;
  occa::memory o_phiRecvBuffer;
  occa::memory o_gPhiSendBuffer;
  occa::memory o_gPhiRecvBuffer;
  occa::memory o_gatherTmpPinned;


  int Nsubsteps;  
  dfloat *Ud, *Ue, *HPhiE, *resU, *rhsUd, sdt;
  occa::memory o_Ud, o_Ue, o_HPhiE, o_resU, o_rhsUd;

  dfloat *cU, *cUd;
  occa::memory o_cU, o_cUd;

  // Some Iso-surfacing variables
  int isoField, isoColorField, isoNfields, isoNlevels, isoMaxNtris, *isoNtris; 
  dfloat isoMinVal, isoMaxVal, *isoLevels, *isoq; 
  size_t isoMax; 
  
  int *isoGNlevels, isoGNgroups;
  dfloat **isoGLvalues;
  // NBN: add storage for compacted isosurf data for gmsh write
  std::vector<dfloat> iso_nodes;
  std::vector<int> iso_tris;


  int readRestartFile,writeRestartFile, restartedFromFile;



  occa::memory *o_isoGLvalues; 
  occa::memory o_isoLevels, o_isoq, o_isoNtris; 
  occa::memory o_plotInterp, o_plotEToV; 



  

  occa::kernel scaledAddKernel;
  occa::kernel multiplyScalarKernel;

  occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
  occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
  occa::kernel subCycleRKUpdateKernel;
  occa::kernel subCycleExtKernel;

  occa::memory o_GSave;

  occa::memory o_U, o_P, o_Phi, o_Psi;
  occa::memory o_Rho, o_Mu, o_GMu;
  occa::memory o_lapPhi;   // laplace of extrapolated phi for  Helmholtz-Psi Rhs  

  occa::memory o_rhsU, o_rhsV, o_rhsW, o_rhsP, o_rhsPhi; 

  occa::memory o_NU, o_LU, o_GP, o_NPhi, o_HPhi, o_GHPhi, o_GPhi;
  occa::memory o_GU, o_DU, o_SU; // o_GU is depreceated

  occa::memory o_UH, o_VH, o_WH, o_PhiH, o_Uhat;
  occa::memory o_rkU, o_rkP, o_rkPhi, o_PI;
  occa::memory o_rkNU, o_rkLU, o_rkGP;

  occa::memory o_Vort, o_Div;

  occa::memory o_vHaloBuffer, o_pHaloBuffer, o_phiHaloBuffer, o_gPhiHaloBuffer; 
  occa::memory o_velocityHaloGatherTmp;

  //ARK data
  occa::memory o_rkC;
  occa::memory o_erkA, o_irkA, o_prkA;
  occa::memory o_erkB, o_irkB, o_prkB;
  occa::memory o_erkE, o_irkE, o_prkE;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

  occa::kernel phaseFieldExtrapolateKernel;  
  occa::kernel phaseFieldAdvectionVolumeKernel;
  occa::kernel phaseFieldAdvectionSurfaceKernel;
  occa::kernel phaseFieldDivGradKernel; 
  occa::kernel phaseFieldUpdateMixingEnergyKernel; 
  occa::kernel phaseFieldAxGradKernel;
  occa::kernel phaseFieldAxGradVolumeKernel;
  occa::kernel phaseFieldAxGradSurfaceKernel;
  occa::kernel phaseFieldAxIpdgKernel;
  occa::kernel phaseFieldAxKernel;


  occa::kernel phaseFieldRhsSolve1Kernel; 
  occa::kernel phaseFieldRhsSolve2Kernel; 
  occa::kernel phaseFieldRhsIpdgBCKernel;
  occa::kernel phaseFieldRhsBCKernel;
  occa::kernel phaseFieldAddBCKernel;

  occa::kernel phaseFieldGradientVolumeKernel;
  occa::kernel phaseFieldGradientSurfaceKernel;


  occa::kernel velocityHaloExtractKernel;
  occa::kernel velocityHaloScatterKernel;
  occa::kernel pressureHaloExtractKernel;
  occa::kernel pressureHaloScatterKernel;
  occa::kernel phaseFieldHaloExtractKernel;
  occa::kernel phaseFieldHaloScatterKernel;
  occa::kernel gradPhaseFieldHaloExtractKernel;
  occa::kernel gradPhaseFieldHaloScatterKernel;

  occa::kernel setFlowFieldKernel;
  occa::kernel setPhaseFieldKernel;
  occa::kernel setMaterialPropertyKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;

  occa::kernel diffusionKernel;
  occa::kernel diffusionIpdgKernel;
  occa::kernel velocityGradientKernel;

  occa::kernel pressureGradientVolumeKernel;
  occa::kernel pressureGradientSurfaceKernel;


  occa::kernel explicitDiffusiveKernel;
  occa::kernel explicitDiffusiveVolumeKernel;
  occa::kernel explicitDiffusiveSurfaceKernel;
  occa::kernel explicitUpdateKernel;

  occa::kernel divergenceVolumeKernel;
  occa::kernel divergenceSurfaceKernel;
  
  
  occa::kernel pressureRhsKernel;
  occa::kernel pressureRhsIpdgBCKernel;
  occa::kernel pressureRhsBCKernel;
  occa::kernel pressureAddBCKernel;
  occa::kernel pressurePenaltyKernel;
  occa::kernel pressureUpdateKernel;

  // Intermediate Velocity Kernels
  occa::kernel velocityExtrapolateKernel;  
  occa::kernel velocityGradientVolumeKernel;
  occa::kernel velocityGradientSurfaceKernel;
  occa::kernel velocityCurlVolumeKernel;
  occa::kernel velocityCurlSurfaceKernel;

  occa::kernel velocityAddPressureKernel;
  occa::kernel velocityAddDiffusiveKernel;
  occa::kernel velocityAddTensionKernel;
  occa::kernel velocityAddStressKernel;
  
  // Velocity Solve Kernels
  occa::kernel velocityRhsKernel;
  occa::kernel velocityRhsIpdgBCKernel;
  occa::kernel velocityRhsBCKernel;
  occa::kernel velocityAddBCKernel;
  occa::kernel velocityUpdateKernel;  
  
  occa::kernel vorticityKernel;
  occa::kernel isoSurfaceKernel;


}mppf_t;

mppf_t *mppfSetup(mesh_t *mesh, setupAide options);

// void insRunARK(ins_t *ins);
// void insRunEXTBDF(ins_t *ins);

 void mppfPlotVTU(mppf_t *mppf, char *fileNameBase);
 void mppfRun(mppf_t *mppf);

 void mppfCahnHilliardRhs(mppf_t *mppf, dfloat time);
 void mppfCahnHilliardSolve(mppf_t *mppf, dfloat time);
 void mppfCahnHilliardUpdate(mppf_t *mppf, dfloat time);


void mppfReport(mppf_t *mppf, dfloat time,  int tstep);
void mppfError(mppf_t *mppf, dfloat time);

void mppfAdvection(mppf_t *mppf, dfloat time);
void mppfPressureRhs  (mppf_t *mppf, dfloat time);
void mppfVelocityRhs  (mppf_t *mppf, dfloat time);

void mppfPressureGradient (mppf_t *mppf, dfloat time);
void mppfVelocityGradient (mppf_t *mppf, dfloat time, occa::memory o_Ue, occa::memory o_GU);

void mppfExplicitDiffusive(mppf_t *mppf, dfloat time, occa::memory o_U, occa::memory o_DU, occa::memory o_SU);



void mppfPressureSolve(mppf_t *mppf, dfloat time, occa::memory o_rkP);
void mppfVelocitySolve(mppf_t *mppf, dfloat time, occa::memory o_rkU);


void mppfDivergence(mppf_t *mppf, dfloat time, occa::memory o_U, occa::memory o_DU);

