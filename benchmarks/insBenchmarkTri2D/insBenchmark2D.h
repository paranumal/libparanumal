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
#include "ellipticTri2D.h"

#define UXID 0
#define UYID 1
#define PRID 2

typedef struct {

  mesh_t *mesh;
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



  dfloat *Ut, *Vt, *Pt, *WN; // For stiffly stable scheme
  dfloat dtfactor;



  dfloat g[2];      // gravitational Acceleration

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

  occa::memory o_vHaloBuffer, o_pHaloBuffer, o_tHaloBuffer;

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

  occa::kernel updateUpdateKernel;

}ins_t;

ins_t *insSetup2D(mesh2D *mesh, int factor, char * options,int Nblocks, int Nnodes,
                  char *vSolverOptions, char *vParAlmondOptions,
                  char *pSolverOptions, char *pParAlmondOptions,
                  char *boundaryHeaderFileName, occa::kernelInfo &kernelInfo);

void insRunBenchmark2D(ins_t *ins, char *options, occa::kernelInfo kernelInfo, char *kernelFileName, int Nblocks, int Nnodes);



