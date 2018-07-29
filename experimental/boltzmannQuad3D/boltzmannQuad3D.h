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


#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh3D.h"

typedef struct {
  
  mesh_t   *mesh;

  char *options;
  
  // INS SOLVER OCCA VARIABLES
  dfloat rho, nu;

  dfloat dt;          // time step
  dfloat lambda;      // helmhotz solver -lap(u) + lamda u
  dfloat finalTime;   // final time to run acoustics to
  int   NtimeSteps;  // number of time steps 
  int   Nstages;     // Number of history states to store
  int   index;       // Index of current state
  int   errorStep; 

  dfloat a0, a1, a2, b0, b1, b2, c0, c1, c2, g0, tau; 

#if 0
  dfloat *resQ;
  dfloat *rhsQ;
  dfloat *Q;
#endif
  
  int Nsubsteps;  

  //  occa::memory o_Q, o_resQ, o_rhsQ;
  
  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel rkUpdateKernel;

  occa::memory o_haloBuffer;

  occa::kernel haloExtractKernel;
  occa::kernel haloScatterKernel;

}solver_t;


solver_t *boltzmannSetupQuad3D(mesh_t *mesh);

void boltzmannRunQuad3D(solver_t *solver);

void boltzmannPlotVTUQuad3D(mesh_t *mesh, char *fileNameBase, int fld);
void boltzmannPlotVTUQuad3DV2(mesh_t *mesh, char *fileNameBase, int tstep);

void boltzmannOccaSetupQuad3D(mesh_t *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);

