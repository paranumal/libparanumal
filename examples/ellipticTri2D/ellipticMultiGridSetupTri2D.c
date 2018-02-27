#include "ellipticTri2D.h"

void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options);
dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);

void AxTri2D(void **args, occa::memory &o_x, occa::memory &o_Ax) {

  solver_t *solver = (solver_t *) args[0];
  dfloat *lambda = (dfloat *) args[1];
  char *options = (char *) args[2];

  ellipticOperator2D(solver,*lambda,o_x,o_Ax,options);
}

void coarsenTri2D(void **args, occa::memory &o_x, occa::memory &o_Rx) {

  solver_t *solver = (solver_t *) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_R = solver->o_R;

  precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);
}

void prolongateTri2D(void **args, occa::memory &o_x, occa::memory &o_Px) {

  solver_t *solver = (solver_t *) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_R = solver->o_R;

  precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
}

void ellipticGather(void **args, occa::memory &o_x, occa::memory &o_Gx) {

  solver_t *solver = (solver_t *) args[0];
  ogs_t *ogs       = (ogs_t *) args[1];
  occa::memory *o_s= (occa::memory *) args[2];
  char *options    = (char *) args[3];
  mesh_t *mesh     = solver->mesh;

  if (strstr(options,"SPARSE")) {
    solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_x, mesh->o_mapSgn, *o_s);
    meshParallelGather(mesh, ogs, *o_s, o_Gx);
  } else {
    meshParallelGather(mesh, ogs, o_x, o_Gx);  
  }
  solver->dotMultiplyKernel(ogs->Ngather, ogs->o_gatherInvDegree, o_Gx, o_Gx);
}

void ellipticScatter(void **args, occa::memory &o_x, occa::memory &o_Sx) {

  solver_t *solver = (solver_t *) args[0];
  ogs_t *ogs       = (ogs_t *) args[1];
  occa::memory *o_s= (occa::memory *) args[2];
  char *options    = (char *) args[3];
  mesh_t *mesh     = solver->mesh;

  
  if (strstr(options,"SPARSE")) {
    meshParallelScatter(mesh, ogs, o_x, *o_s);
    solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, *o_s, mesh->o_mapSgn, o_Sx);
  } else {
    meshParallelScatter(mesh, ogs, o_x, o_Sx);  
  }
}

void buildCoarsenerTri2D(solver_t* solver, mesh2D **meshLevels, int Nf, int Nc, const char* options);

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon,
                                dfloat tau, dfloat lambda, int *BCType,
                                const char *options, const char *parAlmondOptions) {

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;

  //read all the nodes files and load them in a dummy mesh array
  mesh2D **meshLevels = (mesh2D**) calloc(mesh->N+1,sizeof(mesh2D*));
  for (int n=1;n<mesh->N+1;n++) {
    meshLevels[n] = (mesh2D *) calloc(1,sizeof(mesh2D));
    meshLevels[n]->Nverts = mesh->Nverts;
    meshLevels[n]->Nfaces = mesh->Nfaces;
    meshLoadReferenceNodesTri2D(meshLevels[n], n);
  }

  //set the number of MG levels and their degree
  int numLevels;
  int *levelDegree;

  if (strstr(options,"ALLDEGREES")) {
    numLevels = mesh->N;
    levelDegree= (int *) calloc(numLevels,sizeof(int));
    for (int n=0;n<numLevels;n++) levelDegree[n] = mesh->N - n; //all degrees
  } else if (strstr(options,"HALFDEGREES")) {
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
    while (dofs>3) {
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
    while (dofs>3) {
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

  //storage for restriction matrices
  dfloat **R = (dfloat **) calloc(numLevels,sizeof(dfloat*));
  occa::memory *o_R = (occa::memory *) calloc(numLevels,sizeof(occa::memory));

  //initialize parAlmond
  precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);
  agmgLevel **levels = precon->parAlmond->levels;

  //build a solver struct for every degree
  solver_t **solversN = (solver_t**) calloc(mesh->N+1,sizeof(solver_t*));
  solversN[mesh->N] = solver; //top level
  for (int n=1;n<numLevels;n++) {  //build solver for this degree
    int Nf = levelDegree[n-1];
    int Nc = levelDegree[n];
    printf("=============BUIDLING MULTIGRID LEVEL OF DEGREE %d==================\n", Nc);
    solversN[Nc] = ellipticBuildMultigridLevelTri2D(solver,Nc,Nf,BCType,options);
  }

  // set multigrid operators for fine levels
  for (int n=0;n<numLevels-1;n++) {
    int N = levelDegree[n];
    solver_t *solverL = solversN[N];

    //add the level manually
    precon->parAlmond->numLevels++;
    levels[n] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
    levels[n]->gatherLevel = false;   //dont gather this level
    if (strstr(options,"CONTINUOUS")) {//use weighted inner products
      precon->parAlmond->levels[n]->weightedInnerProds = true;
      precon->parAlmond->levels[n]->o_weight = solverL->o_invDegree;
      precon->parAlmond->levels[n]->weight = solverL->invDegree;
    }

    //use the matrix free Ax
    levels[n]->AxArgs = (void **) calloc(3,sizeof(void*));
    levels[n]->AxArgs[0] = (void *) solverL;
    levels[n]->AxArgs[1] = (void *) vlambda;
    levels[n]->AxArgs[2] = (void *) options;
    levels[n]->device_Ax = AxTri2D;

    levels[n]->smoothArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->smoothArgs[0] = (void *) solverL;
    levels[n]->smoothArgs[1] = (void *) levels[n];

    levels[n]->Nrows = mesh->Nelements*solverL->mesh->Np;
    levels[n]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*solverL->mesh->Np;

    if (strstr(options,"CHEBYSHEV")) {
      levels[n]->device_smooth = smoothChebyshevTri2D;

      levels[n]->smootherResidual = (dfloat *) calloc(levels[n]->Ncols,sizeof(dfloat));

      // extra storage for smoothing op
      levels[n]->o_smootherResidual = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
      levels[n]->o_smootherResidual2 = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
      levels[n]->o_smootherUpdate = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
    } else {
      levels[n]->device_smooth = smoothTri2D;

      // extra storage for smoothing op
      levels[n]->o_smootherResidual = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat));
    }

    levels[n]->smootherArgs = (void **) calloc(3,sizeof(void*));
    levels[n]->smootherArgs[0] = (void *) solverL;
    levels[n]->smootherArgs[1] = (void *) options;
    levels[n]->smootherArgs[2] = (void *) vlambda;

    dfloat rateTolerance;    // 0 - accept not approximate patches, 1 - accept all approximate patches
    if(strstr(options, "EXACT")){
      rateTolerance = 0.0;
    } else {
      rateTolerance = 1.0;
    }

    //set up the fine problem smoothing
    if(strstr(options, "OVERLAPPINGPATCH")){
      ellipticSetupSmootherOverlappingPatch(solverL, solverL->precon, levels[n], tau, lambda, BCType, options);
    } else if(strstr(options, "FULLPATCH")){
      ellipticSetupSmootherFullPatch(solverL, solverL->precon, levels[n], tau, lambda, BCType, rateTolerance, options);
    } else if(strstr(options, "FACEPATCH")){
      ellipticSetupSmootherFacePatch(solverL, solverL->precon, levels[n], tau, lambda, BCType, rateTolerance, options);
    } else if(strstr(options, "LOCALPATCH")){
      ellipticSetupSmootherLocalPatch(solverL, solverL->precon, levels[n], tau, lambda, BCType, rateTolerance, options);
    } else { //default to damped jacobi
      ellipticSetupSmootherDampedJacobi(solverL, solverL->precon, levels[n], tau, lambda, BCType, options);
    }
  }

  //report top levels
  if (strstr(options,"VERBOSE")) {
    if((rank==0)&&(numLevels>1)) { //report the upper multigrid levels
      printf("------------------Multigrid Report---------------------\n");
      printf("-------------------------------------------------------\n");
      printf("level|  Degree  |    dimension   |      Smoother       \n");
      printf("     |  Degree  |  (min,max,avg) |      Smoother       \n");
      printf("-------------------------------------------------------\n");
    }

    for(int lev=0; lev<numLevels-1; lev++){

      dlong Nrows = levels[lev]->Nrows;
      hlong hNrows = (hlong) levels[lev]->Nrows;

      dlong minNrows=0, maxNrows=0;
      hlong totalNrows=0;
      dfloat avgNrows;

      MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);
      avgNrows = (dfloat) totalNrows/size;

      if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
      MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, MPI_COMM_WORLD);

      char *smootherString;
      if(strstr(options, "OVERLAPPINGPATCH")){
        smootherString = strdup("OVERLAPPINGPATCH");
      } else if(strstr(options, "FULLPATCH")){
        smootherString = strdup("FULLPATCH");
      } else if(strstr(options, "FACEPATCH")){
        smootherString = strdup("FACEPATCH");
      } else if(strstr(options, "LOCALPATCH")){
        smootherString = strdup("LOCALPATCH");
      } else { //default to damped jacobi
        smootherString = strdup("DAMPEDJACOBI");
      }

      char *smootherOptions1 = strdup(" ");
      char *smootherOptions2 = strdup(" ");
      if (strstr(options,"EXACT")) {
        smootherOptions1 = strdup("EXACT");
      }
      if (strstr(options,"CHEBYSHEV")) {
        smootherOptions2 = strdup("CHEBYSHEV");
      }

      if (rank==0){
        printf(" %3d |   %3d    |    %10.2f  |   %s  \n",
          lev, levelDegree[lev], (dfloat)minNrows, smootherString);
        printf("     |          |    %10.2f  |   %s %s  \n", (dfloat)maxNrows, smootherOptions1, smootherOptions2);
        printf("     |          |    %10.2f  |   \n", avgNrows);
      }
    }
    if((rank==0)&&(numLevels>1)) 
      printf("-------------------------------------------------------\n");
  }

  /* build degree 1 problem and pass to AMG */
  nonZero_t *coarseA;
  long long int nnzCoarseA;
  ogs_t *coarseogs;

  solver_t* solverL = solversN[1];
  int basisNp = solverL->mesh->Np;
  dfloat *basis = NULL;

  if (strstr(options,"BERN")) basis = solverL->mesh->VB;

  hlong *coarseGlobalStarts = (hlong*) calloc(size+1, sizeof(hlong));

  if (strstr(options,"IPDG")) {
    ellipticBuildIpdgTri2D(solverL->mesh, basisNp, basis, tau, lambda, BCType, &coarseA, &nnzCoarseA,coarseGlobalStarts, options);
  } else if (strstr(options,"BRDG")) {
    ellipticBuildBRdgTri2D(solverL->mesh, basisNp, basis, tau, lambda, BCType, &coarseA, &nnzCoarseA,coarseGlobalStarts, options);
  } else if (strstr(options,"CONTINUOUS")) {
    ellipticBuildContinuousTri2D(solverL,lambda,&coarseA,&nnzCoarseA,&coarseogs,coarseGlobalStarts,options);
  }

  hlong *Rows = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  hlong *Cols = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

  for (long long int i=0;i<nnzCoarseA;i++) {
    Rows[i] = coarseA[i].row;
    Cols[i] = coarseA[i].col;
    Vals[i] = coarseA[i].val;
  }

  // build amg starting at level 1
  parAlmondAgmgSetup(precon->parAlmond,
                     coarseGlobalStarts,
                     nnzCoarseA,
                     Rows,
                     Cols,
                     Vals,
                     solver->allNeumann,
                     solver->allNeumannPenalty);
  free(coarseA); free(Rows); free(Cols); free(Vals);

  //tell parAlmond to gather this level
  agmgLevel *coarseLevel = precon->parAlmond->levels[numLevels-1];
  if (strstr(options,"CONTINUOUS")) {
    coarseLevel->gatherLevel = true;
    coarseLevel->weightedInnerProds = false;
    
    coarseLevel->Srhs = (dfloat*) calloc(solverL->mesh->Np*solverL->mesh->Nelements,sizeof(dfloat));
    coarseLevel->Sx   = (dfloat*) calloc(solverL->mesh->Np*solverL->mesh->Nelements,sizeof(dfloat));
    coarseLevel->o_Srhs = solverL->mesh->device.malloc(solverL->mesh->Np*solverL->mesh->Nelements*sizeof(dfloat));
    coarseLevel->o_Sx   = solverL->mesh->device.malloc(solverL->mesh->Np*solverL->mesh->Nelements*sizeof(dfloat));

    coarseLevel->gatherArgs = (void **) calloc(4,sizeof(void*));  
    coarseLevel->gatherArgs[0] = (void *) solverL;
    coarseLevel->gatherArgs[1] = (void *) coarseogs;
    coarseLevel->gatherArgs[2] = (void *) &(coarseLevel->o_Sx);
    coarseLevel->gatherArgs[3] = (void *) options;
    coarseLevel->scatterArgs = coarseLevel->gatherArgs;

    coarseLevel->device_gather  = ellipticGather;
    coarseLevel->device_scatter = ellipticScatter;        
  }

  /* build coarsening and prologation operators to connect levels */
  for(int n=1; n<numLevels; n++) {
    //build coarsen and prologation ops
    int Nf = levelDegree[n-1]; //higher degree
    int Nc = levelDegree[n];  

    solver_t *solverL = solversN[Nc];
    buildCoarsenerTri2D(solverL, meshLevels, Nf, Nc, options);
    
    levels[n]->coarsenArgs = (void **) calloc(1,sizeof(void*));
    levels[n]->coarsenArgs[0] = (void *) solverL;

    levels[n]->prolongateArgs = levels[n]->coarsenArgs;
    
    levels[n]->device_coarsen = coarsenTri2D;
    levels[n]->device_prolongate = prolongateTri2D;
  }

  for (int n=1;n<mesh->N+1;n++) free(meshLevels[n]);
  free(meshLevels);
}



void buildCoarsenerTri2D(solver_t* solver, mesh2D **meshLevels, int Nf, int Nc, const char* options) {

  int NpFine   = meshLevels[Nf]->Np;
  int NpCoarse = meshLevels[Nc]->Np;
  dfloat *P    = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
  dfloat *Ptmp = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));

  //initialize P as identity (which it is for SPARSE)
  for (int i=0;i<NpCoarse;i++) P[i*NpCoarse+i] = 1.0;

  if (strstr(options,"IPDG")) {
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

    if (strstr(options,"BERN")) {
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
  }

  //the coarsen matrix is P^T
  solver->R = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
  for (int i=0;i<NpCoarse;i++) {
    for (int j=0;j<NpFine;j++) {
      solver->R[i*NpFine+j] = P[j*NpCoarse+i];
    }
  }
  solver->o_R = solver->mesh->device.malloc(NpFine*NpCoarse*sizeof(dfloat), solver->R);

  free(P); free(Ptmp);
}