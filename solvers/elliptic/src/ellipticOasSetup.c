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

  /* build one ring patch extension using a single process MPI sub-communicator
     and store in elliptic->precon->ellipticOneRing */
  ellipticBuildOneRing(elliptic, lambda, kernelInfo);

  /* build degree 1 problem and pass to AMG */
  nonZero_t *coarseA;
  dlong nnzCoarseA;
  ogs_t *coarseogs;

  //set up the base level
  elliptic_t* ellipticCoarse;
  if (mesh->N>1) { // assume 
    int Nc = 1;
    int Nf = mesh->N;
    printf("=============BUILDING OAS COARSE LEVEL OF DEGREE %d==================\n", Nc);
    ellipticCoarse = ellipticBuildMultigridLevel(elliptic,Nc,Nf);
  }else{
    ellipticCoarse = elliptic;
  }

  int basisNp = ellipticCoarse->mesh->Np;

  dfloat *basis = NULL;
  
  if (options.compareArgs("BASIS","BERN")) basis = ellipticCoarse->mesh->VB;
  
  hlong *coarseGlobalStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  
  if (options.compareArgs("DISCRETIZATION","IPDG")) {
    ellipticBuildIpdg(ellipticCoarse, basisNp, basis, lambda, &coarseA, &nnzCoarseA,coarseGlobalStarts);
  } else if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ellipticBuildContinuous(ellipticCoarse,lambda,&coarseA,&nnzCoarseA,&coarseogs,coarseGlobalStarts);
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
}
