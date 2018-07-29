/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef BNS_H
#define BNS_H 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h> 

#include "mpi.h"

// #include "mesh.h"
#include "mesh2D.h"
#include "mesh3D.h"

// Block size of reduction 
#define blockSize 256

typedef struct{
  int dim; 
  int elementType;

  int Nfields; // Number of fields

  hlong totalElements; 
  dlong Nblock;
	
  int NrkStages; 
  int frame; 
  int fixed_dt;
	


  int Nrhs;   // Number of RHS 
  int Nimex;  // IMEX order
	
  int Ntscale; // Dummy remove later!!!
  dfloat maxError; 


  dfloat dt, cdt, dtMIN;    // time step
  dfloat startTime ; // Start Time
  dfloat finalTime;  // final time 
  dfloat Lambda2;    // Penalty flux
  dfloat cfl; 

  int NtimeSteps;  // number of time steps
  int Nrk;
  int shiftIndex;    // Rhs index shifting for time steppers
  int fexplicit; 	//Set time stepper type, fully explicit or semi-analytic / imex 
  int pmlcubature;  // Set the cunature integration rule for sigma terms in pml 

  int probeFlag; 
  int errorFlag;
  int reportFlag;
  int writeRestartFile, readRestartFile; 

  int pmlFlag;
  int errorStep;   // number of steps between error calculations
  int reportStep;  // number of steps between error calculations

  int procid; 

  dfloat RT, sqrtRT, tauInv, Ma, Re, nu; // Flow parameters


  mesh_t *mesh; 

  dfloat *q, *rhsq, *resq;
   
  int tmethod; 

  // IMEXRK - Kennedy-Carpanter
  dfloat *rhsqim, *rhsqex, *rkrhsqim, *rkrhsqex; 

  // Pml
  int pmlOrder ; 
  dfloat  sigmaXmax, sigmaYmax, sigmaZmax; 
  dfloat *pmlSigmaX, *pmlSigmaY, *pmlSigmaZ; 
  dfloat *pmlqx, *pmlqy, *pmlqz;
  dfloat *pmlrhsqx, *pmlrhsqy, *pmlrhsqz;
  dfloat *pmlresqx, *pmlresqy, *pmlresqz;

  // Some Iso-surfacing variables
  int isoField, isoColorField, isoNfields, isoNlevels, isoMaxNtris, *isoNtris; 
  dfloat isoMinVal, isoMaxVal, *isoLevels, *isoq; 
  size_t isoMax; 

  occa::memory o_isoLevels, o_isoq, o_isoNtris; 
  occa::memory o_plotInterp, o_plotEToV; 


  // IMEX Coefficients
  dfloat LSIMEX_B[4], LSIMEX_C[4], LSIMEX_ABi[4], LSIMEX_ABe[4], LSIMEX_Ad[4];
  // MRSAAB Coefficients
  dfloat *MRSAAB_A, *MRSAAB_B, *MRSAAB_C, *MRAB_A, *MRAB_B, *MRAB_C;
  // SARK and RK3 Coefficients
  dfloat RK_A[5][5], RK_B[5], RK_C[5], SARK_A[5][5], SARK_B[5], SARK_C[5]; 
  // DOPRI45
  // Change this later this shouldnt be hard-coded
  // dfloat rkC[7], rkA[7*7], rkE[7];
  // dfloat rkC[5], rkA[5*5], rkE[5];

  dfloat *rkC, *rkA, *rkE; 
    
  // IMEXRK
  dfloat *rkCex, *rkAex, *rkBex, *rkEex;
  dfloat *rkCim, *rkAim, *rkBim, *rkEim;



  int *isoGNlevels, isoGNgroups;
  dfloat **isoGLvalues;

  occa::memory *o_isoGLvalues; 




  // NBN: add storage for compacted isosurf data for gmsh write
  std::vector<double> iso_nodes;
  std::vector<int> iso_tris;

  int emethod; 
  int tstep, atstep, rtstep, tstepAccepted, rkp;
  dfloat ATOL, RTOL, time; 
  dfloat *ehist, *dthist; 

  dfloat outputInterval, nextOutputTime;
  int outputForceStep;
  
  int Nvort;     // Number of vorticity fields i.e. 3 or 4 
  dfloat *Vort, *VortMag; 
  occa::memory o_Vort, o_VortMag;


  dfloat *rkq, *rkrhsq, *rkerr, *errtmp;
  dfloat *rkqx, *rkrhsqx;
  dfloat *rkqy, *rkrhsqy;
  dfloat *rkqz, *rkrhsqz;

  occa::memory o_rkq, o_rkrhsq, o_rkerr;
  occa::memory o_errtmp;

  // IMEXRK 
  occa::memory o_rhsqim, o_rhsqex, o_rkrhsqim, o_rkrhsqex;
  occa::memory o_rkAex, o_rkEex, o_rkBex;  
  occa::memory o_rkAim, o_rkEim, o_rkBim; 


  occa::memory o_sendBufferPinned;
  occa::memory o_recvBufferPinned;



  dfloat *fQM; 
  occa::memory o_q,o_rhsq, o_resq, o_fQM;
  occa::memory o_rkA, o_rkE, o_sarkC, o_sarkA, o_sarkE; 
  occa::memory o_rkqx, o_rkqy, o_rkqz, o_rkrhsqx, o_rkrhsqy, o_rkrhsqz; 
  occa::memory o_saveq, o_saveqx, o_saveqy, o_saveqz; // for output minor step of addaptive RK 

  // LS Imex vars
  occa::memory o_qY,   o_qZ,   o_qS;
  occa::memory o_qYx,  o_qZx,  o_qSx;
  occa::memory o_qYy,  o_qZy,  o_qSy;

  occa::memory o_pmlSigmaX, o_pmlSigmaY, o_pmlSigmaZ;
  occa::memory o_pmlqx, o_pmlqy, o_pmlqz; 
  occa::memory o_pmlrhsqx, o_pmlrhsqy, o_pmlrhsqz;
  occa::memory o_pmlresqx, o_pmlresqy, o_pmlresqz;

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel updateKernel;
  occa::kernel traceUpdateKernel;
  occa::kernel relaxationKernel;

  occa::kernel pmlVolumeKernel;
  occa::kernel pmlSurfaceKernel;
  occa::kernel pmlRelaxationKernel;
  occa::kernel pmlUpdateKernel;
  occa::kernel pmlTraceUpdateKernel;

  occa::kernel vorticityKernel;

  occa::kernel isoSurfaceKernel;

  // Boltzmann Imex Kernels
  occa::kernel implicitUpdateKernel;
  occa::kernel pmlImplicitUpdateKernel;
  occa::kernel implicitSolveKernel;
  occa::kernel pmlImplicitSolveKernel;
  //
  occa::kernel residualUpdateKernel;
  occa::kernel pmlResidualUpdateKernel;

  occa::kernel updateStageKernel;
  occa::kernel pmlUpdateStageKernel;

  occa::kernel errorEstimateKernel;

  occa::kernel dotMultiplyKernel; 
        

  // IMEXRK Damping Terms
  occa::kernel pmlDampingKernel; 


}bns_t;


bns_t *bnsSetup(mesh_t *mesh, setupAide &options);

// Pml setup for single rate time discretization
void bnsPmlSetup(bns_t *bns, setupAide &options);

// Pml setup for multi rate time discretization
void bnsMRABPmlSetup(bns_t *bns, setupAide &options);

void bnsRun(bns_t *bns, setupAide &options);
void bnsReport(bns_t *bns, dfloat time, setupAide &options);
void bnsError(bns_t *bns, dfloat time, setupAide &options);
void bnsForces(bns_t *bns, dfloat time, setupAide &options);
void bnsPlotVTU(bns_t *bns, char * FileName);
void bnsIsoPlotVTU(bns_t *bns, int isoNtris, dfloat *isoq, char *fileName);
void bnsIsoWeldPlotVTU(bns_t *bns, char *fileName);

//
void bnsRestartWrite(bns_t *bns, setupAide &options, dfloat time); 
void bnsRestartRead(bns_t *bns, setupAide &options); 
// void bnsRestartSetup(bns_t *bns);


// Function for ramp start
void bnsRampFunction(dfloat t, dfloat *ramp, dfloat *drampdt);

// function for body forcing 
void bnsBodyForce(dfloat t, dfloat *fx, dfloat *fy, dfloat *fz,
		  dfloat *intfx, dfloat *intfy, dfloat *intfz);


// Time stepper coefficients
void bnsTimeStepperCoefficients(bns_t *bns, setupAide &options);
void bnsSAADRKCoefficients(bns_t *bns, setupAide &options);

void bnsMRSAABStep(bns_t *bns, int tstep, int haloBytes,
		   dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options);

void bnsLSERKStep(bns_t *bns, int tstep, int haloBytes,
		  dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options);

void bnsSARKStep(bns_t *bns, dfloat time, int haloBytes,
		 dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options);

void bnsRunEmbedded(bns_t *bns, int haloBytes, dfloat * sendBuffer,
		    dfloat *recvBuffer, setupAide &options);

// Welding Tris

int bnsWeldTriVerts(bns_t *bns, int isoNtris, double *isoq);

void bnsIsoPlotGmsh(bns_t *bns, int isoNtris, char *fname, int tstep, int N_offset,     
  					int E_offset, int plotnum, double plottime,    bool bBinary, int procid);



#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12


#endif

