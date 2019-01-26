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

#include "elliptic.h"

void ellipticOasSetup(elliptic_t *elliptic, dfloat lambda,
		      occa::properties &kernelInfo) {

  mesh_t *mesh = elliptic->mesh;

  setupAide options = elliptic->options;

  /* STAGE 1: build overlapping extended partition problem */
  
  /* build one ring patch extension using a single process MPI sub-communicator
     and store in elliptic->precon->ellipticOneRing */
  ellipticBuildOneRing(elliptic, lambda, kernelInfo);
  
  /* STAGE 2: build coarse problem */
  nonZero_t *coarseA;
  dlong nnzCoarseA;
  ogs_t *coarseogs;

  //set up the base level
  int Nc = 1;
  int Nf = mesh->N;

  // build coarsener
  int NqFine   = Nf+1;
  int NqCoarse = Nc+1;

  int NpFine   = (Nf+1)*(Nf+1)*(Nf+1);
  int NpCoarse = (Nc+1)*(Nc+1)*(Nc+1);

  int NblockVFine = maxNthreads/NpFine;
  int NblockVCoarse = maxNthreads/NpCoarse;
  
  elliptic_t* ellipticOasCoarse;
  if (mesh->N>1) { // assume 
    printf("=============BUILDING OAS COARSE LEVEL OF DEGREE %d==================\n", Nc);
    ellipticOasCoarse = ellipticBuildMultigridLevel(elliptic,Nc,Nf);
    
  }else{
    ellipticOasCoarse = elliptic;
  }

  
  dfloat *P    = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  dfloat *R    = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));

  // hard wire for linears
  for(int n=0;n<NqFine;++n){
    P[n*NqCoarse + 0] = 0.5*(1-mesh->gllz[n]);
    P[n*NqCoarse + 1] = 0.5*(1+mesh->gllz[n]);
    R[0*NqFine + n] = 0.5*(1-mesh->gllz[n]);
    R[1*NqFine + n] = 0.5*(1+mesh->gllz[n]);
  }
  
  occa::memory o_R = elliptic->mesh->device.malloc(NqFine*NqCoarse*sizeof(dfloat), R);
  occa::memory o_P = elliptic->mesh->device.malloc(NqFine*NqCoarse*sizeof(dfloat), P);

  free(P); free(R);

  int basisNp = ellipticOasCoarse->mesh->Np;

  hlong *coarseGlobalStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ellipticBuildContinuous(ellipticOasCoarse,lambda,&coarseA,&nnzCoarseA,&coarseogs,coarseGlobalStarts);
  }

  hlong *Rows = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  hlong *Cols = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));
  
  for (dlong i=0;i<nnzCoarseA;i++) {
    Rows[i] = coarseA[i].row;
    Cols[i] = coarseA[i].col;
    Vals[i] = coarseA[i].val;
  }

  printf("nnzCoarseA = %d\n", nnzCoarseA);
  
  free(coarseA);

  elliptic->precon->ellipticOasCoarse = ellipticOasCoarse;
  elliptic->precon->o_oasRestrictionMatrix = o_R;
  elliptic->precon->o_oasProlongationMatrix = o_P;
  elliptic->precon->o_oasCoarseTmp = mesh->device.malloc(NpCoarse*mesh->Nelements*sizeof(dfloat));
  elliptic->precon->o_oasFineTmp   = mesh->device.malloc(NpFine*mesh->Nelements*sizeof(dfloat));
  
  // build degree 1 coarsening and prolongation matrices and kernels
  
  kernelInfo["defines/" "p_NqFine"]= Nf+1;
  kernelInfo["defines/" "p_NqCoarse"]= Nc+1;

  kernelInfo["defines/" "p_NpFine"]= NpFine;
  kernelInfo["defines/" "p_NpCoarse"]= NpCoarse;
  
  kernelInfo["defines/" "p_NblockVFine"]= NblockVFine;
  kernelInfo["defines/" "p_NblockVCoarse"]= NblockVCoarse;

  char *suffix;

  if(elliptic->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];
  
  sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
  sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
  elliptic->precon->oasRestrictionKernel =
    mesh->device.buildKernel(fileName,kernelName,kernelInfo);
  
  sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
  sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
  elliptic->precon->oasProlongationKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
  
  
}
