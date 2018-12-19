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

  mesh_t     *mesh;
  elliptic_t *solver;
  int NVfields;            // Number of velocity fields
  int NSfields;            // Number of scalar fields
  
  setupAide options;
  // INS SOLVER OCCA VARIABLES
  dfloat k, nu, cp, rho, ialf, alf,  Re, Pr;
  dfloat ubar, vbar, wbar, sbar;
  dlong vOffset;
  dlong sOffset;
  dlong Ntotal;
  int Nblock;
  dfloat dt, idt, cfl, dti;          // time step
  dfloat dtMIN;         
  dfloat time;
  int tstep, frame;
  dfloat g0, ig0, lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat startTime;   
  dfloat finalTime;   

  int temporalOrder;
  int ExplicitOrder; 
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     
  int   outputStep;
  int   outputForceStep; 
  int   dtAdaptStep; 
  
  int Niter;


  //solver tolerances
  dfloat TOL;

  dfloat icp, irho; // hold some inverses
  
  dfloat *U, *S;
  dfloat *NS, *rkNS;
  //  dfloat *rhsS;   
  dfloat *rkS; 
  dfloat g[3];      // gravitational Acceleration

  //RK Subcycle Data
  int SNrk;
  dfloat *Srka, *Srkb, *Srkc; 
  //EXTBDF data
  dfloat *extbdfA, *extbdfB, *extbdfC;
  dfloat *extC;

  int *mapB;
  occa::memory o_mapB;

  //halo data
  dfloat *sendBuffer;
  dfloat *recvBuffer;
  dfloat *haloGatherTmp;

  occa::memory o_sendBuffer;
  occa::memory o_recvBuffer;
  occa::memory o_gatherTmpPinned;

  int Nsubsteps;
  dfloat sdt; 
  dfloat *Sd, *Ue, *resS, *rhsS, *rhsSd;
  occa::memory o_Sd, o_Ue, o_resS, o_rhsS, o_rhsSd;

  dfloat *cU, *cSd, *cS; 
  occa::memory o_cU, o_cSd, o_cS;

  // Some Iso-surfacing variables
  // int isoField, isoColorField, isoNfields, isoNlevels, isoMaxBNtris, *isoNtris; 
  // dfloat isoMinVal, isoMaxVal, *isoLevels, *isoq; 
  // size_t isoMax; 
  
  // int *isoGNlevels, isoGNgroups;
  // dfloat **isoGLvalues;
  // NBN: add storage for compacted isosurf data for gmsh write
  // std::vector<dfloat> iso_nodes;
  // std::vector<int> iso_tris;
  // int readRestartFile,writeRestartFile, restartedFromFile;
  // occa::memory *o_isoGLvalues; 
  // occa::memory o_isoLevels, o_isoq, o_isoNtris; 
  // occa::memory o_plotInterp, o_plotEToV; 

   occa::kernel scaledAddKernel;
   occa::kernel subCycleVolumeKernel,  subCycleCubatureVolumeKernel ;
   occa::kernel subCycleSurfaceKernel, subCycleCubatureSurfaceKernel;;
   occa::kernel subCycleRKUpdateKernel;
   occa::kernel subCycleExtKernel;

  // occa::kernel constrainKernel;
  
  occa::memory o_U; 
  occa::memory o_S, o_SH, o_NS;
  occa::memory o_rkS, o_rkNS; 
  
  // occa::memory o_Vort, o_Div; // Not sure to keep it
  occa::memory o_haloBuffer;
  occa::memory o_haloGatherTmp;

  //ARK data
  occa::memory o_rkC;
  occa::memory o_erkA, o_irkA, o_prkA;
  occa::memory o_erkB, o_irkB, o_prkB;
  occa::memory o_erkE, o_irkE, o_prkE;

  //EXTBDF data
  occa::memory o_extbdfA, o_extbdfB, o_extbdfC;
  occa::memory o_extC;

  occa::kernel haloExtractKernel;
  occa::kernel haloScatterKernel;

  occa::kernel setFlowFieldKernel;
  occa::kernel setScalarFieldKernel;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;
  occa::kernel advectionCubatureVolumeKernel;
  occa::kernel advectionCubatureSurfaceKernel;

  occa::kernel helmholtzRhsKernel;
  occa::kernel helmholtzRhsIpdgBCKernel;
  occa::kernel helmholtzRhsBCKernel;
  occa::kernel helmholtzAddBCKernel;
    
}cds_t;

cds_t *cdsSetup(mesh_t *mesh, setupAide options);

void cdsRunEXTBDF(cds_t *cds);

void cdsPlotVTU(cds_t *cds, char *fileNameBase);
void cdsReport(cds_t *cds, dfloat time,  int tstep);
void cdsError(cds_t *cds, dfloat time);

void cdsAdvection(cds_t *cds, dfloat time, occa::memory o_U, occa::memory o_S, occa::memory o_NS);
void cdsSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_U, occa::memory o_S, occa::memory o_NU);

void cdsHelmholtzRhs(cds_t *cds, dfloat time, int stage, occa::memory o_rhsS);
void cdsHelmholtzSolve(cds_t *cds, dfloat time, int stage, occa::memory o_rhsS,occa::memory o_rkS);


// // Restarting from file
// void insRestartWrite(ins_t *ins, setupAide &options, dfloat time); 
// void insRestartRead(ins_t *ins, setupAide &options); 

// customized hex writer
// extern "C"
// {
//   void insPlotVTUHex3D(ins_t *ins, char *fileNameBase);
//   void insPlotWallsVTUHex3D(ins_t *ins, char *fileNameBase);
// }


// void insRenderQuad3D(ins_t *ins, char *fileBaseName, int fileIndex);

// void simpleRayTracer(int     plotNelements,
// 		     dfloat *plotx,
// 		     dfloat *ploty,
// 		     dfloat *plotz,
// 		     dfloat *plotq,
// 		     const char *fileBaseName,
// 		     const int fileIndex);
