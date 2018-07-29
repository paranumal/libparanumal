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

void ellipticMultigridAx(void **args, occa::memory &o_x, occa::memory &o_Ax) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  dfloat *lambda = (dfloat *) args[1];

  ellipticOperator(elliptic,*lambda,o_x,o_Ax, dfloatString); // "float" ); // hard coded for testing (should make an option)
}

void ellipticMultigridCoarsen(void **args, occa::memory &o_x, occa::memory &o_Rx) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  elliptic_t *Felliptic = (elliptic_t *) args[1];
  setupAide options = elliptic->options;

  mesh_t *mesh = elliptic->mesh;
  mesh_t *Fmesh = Felliptic->mesh;
  precon_t *precon = elliptic->precon;
  occa::memory o_R = elliptic->o_R;

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS"))
    Felliptic->dotMultiplyKernel(Fmesh->Nelements*Fmesh->Np, Fmesh->ogs->o_invDegree, o_x, o_x);

  precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ellipticParallelGatherScatter(mesh, mesh->ogs, o_Rx, dfloatString, "add");  
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_Rx);
  }
}

void ellipticMultigridProlongate(void **args, occa::memory &o_x, occa::memory &o_Px) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  mesh_t *mesh = elliptic->mesh;
  precon_t *precon = elliptic->precon;
  occa::memory o_R = elliptic->o_R;

  precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
}

void ellipticGather(void **args, occa::memory &o_x, occa::memory &o_Gx) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  ogs_t *ogs       = (ogs_t *) args[1];
  occa::memory *o_s= (occa::memory *) args[2];
  
  mesh_t *mesh      = elliptic->mesh;
  setupAide options = elliptic->options;

  meshParallelGather(mesh, ogs, o_x, o_Gx);  
  elliptic->dotMultiplyKernel(ogs->Ngather, ogs->o_gatherInvDegree, o_Gx, o_Gx);
}

void ellipticScatter(void **args, occa::memory &o_x, occa::memory &o_Sx) {

  elliptic_t *elliptic = (elliptic_t *) args[0];
  ogs_t *ogs       = (ogs_t *) args[1];
  occa::memory *o_s= (occa::memory *) args[2];
  
  mesh_t *mesh      = elliptic->mesh;
  setupAide options = elliptic->options;

  meshParallelScatter(mesh, ogs, o_x, o_Sx);  
}

void buildCoarsenerTriTet(elliptic_t* elliptic, mesh_t **meshLevels, int Nf, int Nc);
void buildCoarsenerQuadHex(elliptic_t* elliptic, mesh_t **meshLevels, int Nf, int Nc);

void ellipticMultiGridSetup(elliptic_t *elliptic, precon_t* precon, dfloat lambda) {

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  //read all the nodes files and load them in a dummy mesh array
  mesh_t **meshLevels = (mesh_t**) calloc(mesh->N+1,sizeof(mesh_t*));
  for (int n=1;n<mesh->N+1;n++) {
    meshLevels[n] = (mesh_t *) calloc(1,sizeof(mesh_t));
    meshLevels[n]->Nverts = mesh->Nverts;
    meshLevels[n]->Nfaces = mesh->Nfaces;
    
    switch(elliptic->elementType){
    case TRIANGLES:
      meshLoadReferenceNodesTri2D(meshLevels[n], n); break;
    case QUADRILATERALS:
      meshLoadReferenceNodesQuad2D(meshLevels[n], n); break;
    case TETRAHEDRA:
      meshLoadReferenceNodesTet3D(meshLevels[n], n); break;
    case HEXAHEDRA:
      meshLoadReferenceNodesHex3D(meshLevels[n], n); break;
    }
  }

  //set the number of MG levels and their degree
  int numLevels;
  int *levelDegree;

  if (options.compareArgs("MULTIGRID COARSENING","ALLDEGREES")) {
    numLevels = mesh->N;
    levelDegree= (int *) calloc(numLevels,sizeof(int));
    for (int n=0;n<numLevels;n++) levelDegree[n] = mesh->N - n; //all degrees
  } else if (options.compareArgs("MULTIGRID COARSENING","HALFDEGREES")) {
    numLevels = floor(mesh->N/2.)+1;
    levelDegree= (int *) calloc(numLevels,sizeof(int));
    for (int n=0;n<numLevels;n++) levelDegree[n] = mesh->N - 2*n; //decrease by two
    levelDegree[numLevels-1] = 1; //ensure the last level is degree 1
  } else { //default "HALFDOFS"
    // pick the degrees so the dofs of each level halfs (roughly)
    //start by counting the number of levels neccessary
    numLevels = 1;
    int degree = mesh->N;
    int dofs = meshLevels[degree]->Np;
    int basedofs = mesh->Nverts;
    while (dofs>basedofs) {
      numLevels++;
      for (;degree>0;degree--)
        if (meshLevels[degree]->Np<=dofs/2)
          break;
      dofs = meshLevels[degree]->Np;
    }
    levelDegree= (int *) calloc(numLevels,sizeof(int));
    degree = mesh->N;
    numLevels = 1;
    levelDegree[0] = degree;
    dofs = meshLevels[degree]->Np;
    while (dofs>basedofs) {
      for (;degree>0;degree--)
        if (meshLevels[degree]->Np<=dofs/2)
          break;
      dofs = meshLevels[degree]->Np;
      levelDegree[numLevels] = degree;
      numLevels++;
    }
  }

  //storage for lambda parameter
  dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
  *vlambda = lambda;

  //initialize parAlmond
  precon->parAlmond = parAlmondInit(mesh, elliptic->options);
  agmgLevel **levels = precon->parAlmond->levels;

  //build a elliptic struct for every degree
  elliptic_t **ellipticsN = (elliptic_t**) calloc(mesh->N+1,sizeof(elliptic_t*));
  ellipticsN[mesh->N] = elliptic; //top level
  for (int n=1;n<numLevels;n++) {  //build elliptic for this degree
    int Nf = levelDegree[n-1];
    int Nc = levelDegree[n];
    printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nc);
    ellipticsN[Nc] = ellipticBuildMultigridLevel(elliptic,Nc,Nf);
  }

  // set multigrid operators for fine levels
  for (int n=0;n<numLevels-1;n++) {
    int N = levelDegree[n];
    elliptic_t *ellipticL = ellipticsN[N];

    //add the level manually
    precon->parAlmond->numLevels++;
    levels[n] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
    levels[n]->gatherLevel = false;   //dont gather this level
    if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {//use weighted inner products
      precon->parAlmond->levels[n]->weightedInnerProds = true;
      precon->parAlmond->levels[n]->o_weight = ellipticL->o_invDegree;
      precon->parAlmond->levels[n]->weight = ellipticL->invDegree;
    }

    //use the matrix free Ax
    levels[n]->AxArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->AxArgs[0] = (void *) ellipticL;
    levels[n]->AxArgs[1] = (void *) vlambda;
    levels[n]->device_Ax = ellipticMultigridAx;

    levels[n]->smoothArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->smoothArgs[0] = (void *) ellipticL;
    levels[n]->smoothArgs[1] = (void *) levels[n];

    levels[n]->Nrows = mesh->Nelements*ellipticL->mesh->Np;
    levels[n]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*ellipticL->mesh->Np;

    if (options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {
      if (!options.getArgs("MULTIGRID CHEBYSHEV DEGREE", levels[n]->ChebyshevIterations))
        levels[n]->ChebyshevIterations = 2; //default to degree 2

      levels[n]->device_smooth = ellipticMultigridSmoothChebyshev;

      levels[n]->smootherResidual = (dfloat *) calloc(levels[n]->Ncols,sizeof(dfloat));

      // extra storage for smoothing op
      levels[n]->o_smootherResidual = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
      levels[n]->o_smootherResidual2 = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
      levels[n]->o_smootherUpdate = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
    } else {
      levels[n]->device_smooth = ellipticMultigridSmooth;

      // extra storage for smoothing op
      levels[n]->o_smootherResidual = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat));
    }

    levels[n]->smootherArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->smootherArgs[0] = (void *) ellipticL;
    levels[n]->smootherArgs[1] = (void *) vlambda;

    dfloat rateTolerance;    // 0 - accept not approximate patches, 1 - accept all approximate patches
    if(options.compareArgs("MULTIGRID SMOOTHER","EXACT")){
      rateTolerance = 0.0;
    } else {
      rateTolerance = 1.0;
    }

    //set up the fine problem smoothing
    if(options.compareArgs("MULTIGRID SMOOTHER","LOCALPATCH")){
      ellipticSetupSmootherLocalPatch(ellipticL, ellipticL->precon, levels[n], lambda, rateTolerance);
    } else { //default to damped jacobi
      ellipticSetupSmootherDampedJacobi(ellipticL, ellipticL->precon, levels[n], lambda);
    }
  }

  //report top levels
  if (options.compareArgs("VERBOSE","TRUE")) {
    if((mesh->rank==0)&&(numLevels>0)) { //report the upper multigrid levels
      printf("------------------Multigrid Report---------------------\n");
      printf("-------------------------------------------------------\n");
      printf("level|  Degree  |    dimension   |      Smoother       \n");
      printf("     |  Degree  |  (min,max,avg) |      Smoother       \n");
      printf("-------------------------------------------------------\n");
    }

    for(int lev=0; lev<numLevels; lev++){

      dlong Nrows = (lev==numLevels-1) ? mesh->Nverts*mesh->Nelements: levels[lev]->Nrows;
      hlong hNrows = (hlong) Nrows;

      dlong minNrows=0, maxNrows=0;
      hlong totalNrows=0;
      dfloat avgNrows;

      MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, mesh->comm);
      MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, mesh->comm);
      avgNrows = (dfloat) totalNrows/mesh->size;

      if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
      MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, mesh->comm);

      char smootherString[BUFSIZ];
      strcpy(smootherString, (char*) (options.getArgs("MULTIGRID SMOOTHER")).c_str());

      if (mesh->rank==0){
        printf(" %3d |   %3d    |    %10.2f  |   %s  \n",
          lev, levelDegree[lev], (dfloat)minNrows, smootherString);
        printf("     |          |    %10.2f  |   \n", (dfloat)maxNrows);
        printf("     |          |    %10.2f  |   \n", avgNrows);
      }
    }
    if((mesh->rank==0)&&(numLevels>0)) 
      printf("-------------------------------------------------------\n");
  }

  /* build degree 1 problem and pass to AMG */
  nonZero_t *coarseA;
  dlong nnzCoarseA;
  ogs_t *coarseogs;

  elliptic_t* ellipticL = ellipticsN[1];
  int basisNp = ellipticL->mesh->Np;
  dfloat *basis = NULL;

  if (options.compareArgs("BASIS","BERN")) basis = ellipticL->mesh->VB;

  hlong *coarseGlobalStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));

  if (options.compareArgs("DISCRETIZATION","IPDG")) {
    ellipticBuildIpdg(ellipticL, basisNp, basis, lambda, &coarseA, &nnzCoarseA,coarseGlobalStarts);
  } else if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ellipticBuildContinuous(ellipticL,lambda,&coarseA,&nnzCoarseA,&coarseogs,coarseGlobalStarts);
  }

  hlong *Rows = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  hlong *Cols = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

  for (dlong i=0;i<nnzCoarseA;i++) {
    Rows[i] = coarseA[i].row;
    Cols[i] = coarseA[i].col;
    Vals[i] = coarseA[i].val;
  }

  // build amg starting at level N=1
  parAlmondAgmgSetup(precon->parAlmond,
                     coarseGlobalStarts,
                     nnzCoarseA,
                     Rows,
                     Cols,
                     Vals,
                     elliptic->allNeumann,
                     elliptic->allNeumannPenalty);
  free(coarseA); free(Rows); free(Cols); free(Vals);

  //tell parAlmond to gather this level
  agmgLevel *coarseLevel = precon->parAlmond->levels[numLevels-1];
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    coarseLevel->gatherLevel = true;
    coarseLevel->weightedInnerProds = false;
    
    coarseLevel->Srhs = (dfloat*) calloc(ellipticL->mesh->Np*ellipticL->mesh->Nelements,sizeof(dfloat));
    coarseLevel->Sx   = (dfloat*) calloc(ellipticL->mesh->Np*ellipticL->mesh->Nelements,sizeof(dfloat));
    coarseLevel->o_Srhs = ellipticL->mesh->device.malloc(ellipticL->mesh->Np*ellipticL->mesh->Nelements*sizeof(dfloat),coarseLevel->Srhs);
    coarseLevel->o_Sx   = ellipticL->mesh->device.malloc(ellipticL->mesh->Np*ellipticL->mesh->Nelements*sizeof(dfloat),coarseLevel->Sx);

    coarseLevel->gatherArgs = (void **) calloc(3,sizeof(void*));  
    coarseLevel->gatherArgs[0] = (void *) ellipticL;
    coarseLevel->gatherArgs[1] = (void *) coarseogs;
    coarseLevel->gatherArgs[2] = (void *) &(coarseLevel->o_Sx);
    coarseLevel->scatterArgs = coarseLevel->gatherArgs;

    coarseLevel->device_gather  = ellipticGather;
    coarseLevel->device_scatter = ellipticScatter;        
  }

  /* build coarsening and prologation operators to connect levels */
  for(int n=1; n<numLevels; n++) {
    //build coarsen and prologation ops
    int Nf = levelDegree[n-1]; //higher degree
    int Nc = levelDegree[n];  

    elliptic_t *ellipticL = ellipticsN[Nc];
    elliptic_t *ellipticF = ellipticsN[Nf];

    if (elliptic->elementType==TRIANGLES||elliptic->elementType==TETRAHEDRA){
      buildCoarsenerTriTet(ellipticL, meshLevels, Nf, Nc);
    } else {
      buildCoarsenerQuadHex(ellipticL, meshLevels, Nf, Nc);
    }
    
    levels[n]->coarsenArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->coarsenArgs[0] = (void *) ellipticL;
    levels[n]->coarsenArgs[1] = (void *) ellipticF;

    levels[n]->prolongateArgs = levels[n]->coarsenArgs;
    
    levels[n]->device_coarsen = ellipticMultigridCoarsen;
    levels[n]->device_prolongate = ellipticMultigridProlongate;
  }

  for (int n=1;n<mesh->N+1;n++) free(meshLevels[n]);
  free(meshLevels);
}



void buildCoarsenerTriTet(elliptic_t* elliptic, mesh_t **meshLevels, int Nf, int Nc) {

  int NpFine   = meshLevels[Nf]->Np;
  int NpCoarse = meshLevels[Nc]->Np;
  dfloat *P    = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
  dfloat *Ptmp = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));

  //initialize P as identity (which it is for SPARSE)
  for (int i=0;i<NpCoarse;i++) P[i*NpCoarse+i] = 1.0;


  for (int n=Nc;n<Nf;n++) {

    int Npp1 = meshLevels[n+1]->Np;
    int Np   = meshLevels[n  ]->Np;

    //copy P
    for (int i=0;i<Np*NpCoarse;i++) Ptmp[i] = P[i];

    //Multiply by the raise op
    for (int i=0;i<Npp1;i++) {
      for (int j=0;j<NpCoarse;j++) {
        P[i*NpCoarse + j] = 0.;
        for (int k=0;k<Np;k++) {
          P[i*NpCoarse + j] += meshLevels[n]->interpRaise[i*Np+k]*Ptmp[k*NpCoarse + j];
        }
      }
    }
  }

  if (elliptic->options.compareArgs("BASIS","BERN")) {
    dfloat* BBP = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
    for (int j=0;j<NpFine;j++) {
      for (int i=0;i<NpCoarse;i++) {
        for (int k=0;k<NpCoarse;k++) {
          for (int l=0;l<NpFine;l++) {
            BBP[i+j*NpCoarse] += meshLevels[Nf]->invVB[l+j*NpFine]*P[k+l*NpCoarse]*meshLevels[Nc]->VB[i+k*NpCoarse];
          }
        }
      }
    }
    for (int j=0;j<NpFine;j++) {
      for (int i=0;i<NpCoarse;i++) {
        P[i+j*NpCoarse] = BBP[i+j*NpCoarse];
      }
    }
    free(BBP);
  }

  //the coarsen matrix is P^T
  elliptic->R = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
  for (int i=0;i<NpCoarse;i++) {
    for (int j=0;j<NpFine;j++) {
      elliptic->R[i*NpFine+j] = P[j*NpCoarse+i];
    }
  }
  elliptic->o_R = elliptic->mesh->device.malloc(NpFine*NpCoarse*sizeof(dfloat), elliptic->R);

  free(P); free(Ptmp);
}

void buildCoarsenerQuadHex(elliptic_t* elliptic, mesh_t **meshLevels, int Nf, int Nc) {

  int NqFine   = Nf+1;
  int NqCoarse = Nc+1;
  dfloat *P    = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  dfloat *Ptmp = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));

  //initialize P as identity
  for (int i=0;i<NqCoarse;i++) P[i*NqCoarse+i] = 1.0;

  for (int n=Nc;n<Nf;n++) {

    int Nqp1 = n+2;
    int Nq   = n+1;

    //copy P
    for (int i=0;i<Nq*NqCoarse;i++) Ptmp[i] = P[i];

    //Multiply by the raise op
    for (int i=0;i<Nqp1;i++) {
      for (int j=0;j<NqCoarse;j++) {
        P[i*NqCoarse + j] = 0.;
        for (int k=0;k<Nq;k++) {
          P[i*NqCoarse + j] += meshLevels[n]->interpRaise[i*Nq+k]*Ptmp[k*NqCoarse + j];
        }
      }
    }
  }

  //the coarsen matrix is P^T
  elliptic->R = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  for (int i=0;i<NqCoarse;i++) {
    for (int j=0;j<NqFine;j++) {
      elliptic->R[i*NqFine+j] = P[j*NqCoarse+i];
    }
  }
  elliptic->o_R = elliptic->mesh->device.malloc(NqFine*NqCoarse*sizeof(dfloat), elliptic->R);

  free(P); free(Ptmp);
}
