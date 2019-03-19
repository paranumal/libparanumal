#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh3D.h"

typedef struct {
  
    mesh_t   *mesh;

    //physics parameters
    iint Nfields;
    dfloat RT;
    dfloat sqrtRT;
    dfloat rho, nu, tau, tauInv;

    //occa common (for now)
    occa::device device;
    char deviceConfig[BUFSIZ];
    occa::memory o_q;
    occa::memory o_rhsq;
    occa::memory o_rhsqs;
    occa::memory o_rhsqw;
    occa::memory o_LIFT;
    occa::memory o_LIFTT;
    occa::memory o_vgeo;
    occa::memory o_sgeo;
    occa::memory o_vmapM;
    occa::memory o_vmapP;
    occa::memory o_mapP;
    occa::memory o_x;
    occa::memory o_y;
    occa::memory o_z;
    occa::memory o_haloElementList;
    occa::memory o_haloBuffer;
    occa::memory o_D;
    occa::memory o_weakD;
    occa::memory o_dualProjMatrix;
    occa::memory o_cubeFaceNumber;
    occa::memory o_EToE;
    occa::memory o_qCorr;
    occa::memory o_qFilter;
    occa::memory o_resq;
    occa::memory o_eInterp;
    occa::memory o_overlapDirection;
    occa::memory o_cubeDistance;
    occa::memory o_rlocal;
    occa::memory o_slocal;
    occa::memory o_par_loc;
    occa::memory o_perp_index;
    occa::memory o_gridToE;
    occa::memory o_mass;
    occa::memory o_invmass;
    occa::kernel haloExtractKernel;
    occa::kernel haloScatterKernel;
    occa::kernel volumeKernel;
    occa::kernel volumeCorrectionKernel;
    occa::kernel surfaceKernel;
    occa::kernel loadFilterGridKernel;
    occa::kernel updateKernel;
    occa::kernel filterKernelH;
    occa::kernel filterKernelV;
    occa::kernel massMatrixKernel;
    
    //occa mrsaab
    occa::memory o_qpre;
    occa::memory o_qPreCorr;
    occa::memory o_prerhsq;
    occa::memory o_qPreFilter;
    occa::memory o_qPreFiltered;
    occa::memory o_MRABlevels;
    occa::memory *o_MRABelementIds;
    occa::memory *o_MRABhaloIds;
    occa::memory o_fQ;
    occa::memory o_lev_updates;
    occa::memory o_shift;
    occa::memory o_qFiltered;
    occa::kernel volumePreKernel;
    occa::kernel volumeCorrPreKernel;
    occa::kernel surfacePreKernel;
    occa::kernel traceUpdateKernel;
    occa::kernel traceUpdatePreKernel;
    occa::kernel updatePreKernel;
    occa::kernel filterKernelq0H;
    occa::kernel filterKernelq0V;
    occa::kernel filterKernelHLSERK;
    occa::kernel filterKernelVLSERK;
    occa::kernel filterKernelLevelsH;
    occa::kernel filterKernelLevelsV;

    //occa dopri
    occa::memory o_rkA;
    occa::memory o_rkE;
    occa::memory o_rkq;
    occa::memory o_rkerr;
    occa::memory o_errtmp;
    occa::kernel rkStageKernel;
    occa::kernel rkErrorEstimateKernel;
  
    //common parameters (for now)
    dfloat *q;
    iint Nrk;
    iint NtimeSteps;
    dfloat *resq;
    dfloat dt;          // time step
    dfloat finalTime;   // final time to run acoustics to
  
    //MRSAAB parameters
    iint Nrhs;
    dfloat *fQ;
    dfloat *rhsq;
    iint *MRABshiftIndex;
    iint *lev_updates;
    dfloat *localdt;
    iint *shift;
    iint errorStep;
    dfloat *MRSAAB_A;
    dfloat *MRSAAB_B;
    dfloat *MRSAAB_C;
    dfloat *MRAB_A;
    dfloat *MRAB_B;
    dfloat *MRAB_C;
    dfloat *rka;
    dfloat *rkb;
    dfloat *rkc;
  
    //DOPRI parameters
    dfloat *rkA;
    dfloat *rkE;
    dfloat *rkC;
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

void advectionRunMRSAABQuad3D(solver_t *solver);
void advectionRunLSERKQuad3D(solver_t *solver);
void advectionRunDOPRIQuad3D(solver_t *solver);
void advectionRunLSERKbasicQuad3D(solver_t *solver,dfloat alpha_scale);
void advectionRunLSERKsymQuad3D(solver_t *solver,dfloat alpha_scale);

void advectionSpectrumLSERKQuad3D(solver_t *solver,dfloat alpha_scale);

void advectionPlotVTUQuad3D(mesh_t *mesh, char *fileNameBase, iint fld);
void advectionPlotVTUQuad3DV2(solver_t *solver, char *fileNameBase, iint tstep);

void advectionSetupOccaQuad3D(solver_t *solver, occa::kernelInfo *kernelInfo);
solver_t *advectionSetupPhysicsQuad3D(mesh_t *mesh);
void advectionSetupMRSAABQuad3D(solver_t *solver);
void advectionSetupDOPRIQuad3D(solver_t *solver);
void advectionSetupLSERKQuad3D(solver_t *solver);
void advectionSetupLSERKsymQuad3D(solver_t *solver);

void meshMRABSetupQuad3D(mesh3D *mesh, dfloat *EToDT, int maxLevels);

//void advectionPlotLevels(mesh_t *mesh, char *fileNameBase, iint tstep,dfloat *q);

//void advectionPlotNorms(mesh_t *mesh, char *fileNameBase, iint tstep,dfloat *q);

void advectionErrorNormQuad3D(solver_t *solver, dfloat t, char *fileBase, int slice);
