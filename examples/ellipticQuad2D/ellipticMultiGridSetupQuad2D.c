#include "ellipticQuad2D.h"

void AxQuad2D(void **args, occa::memory &o_x, occa::memory &o_Ax) {

  solver_t *solver = (solver_t *) args[0];
  dfloat *lambda = (dfloat *) args[1];
  char *options = (char *) args[2];

  ellipticOperator2D(solver,*lambda,o_x,o_Ax,options);
}

void coarsenQuad2D(void **args, occa::memory &o_x, occa::memory &o_Rx) {

  solver_t *solver = (solver_t *) args[0];
  char *options    = (char *) args[1];
  solver_t *Fsolver = (solver_t *) args[2];

  mesh_t *mesh = solver->mesh;
  mesh_t *Fmesh = Fsolver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_R = solver->o_R;

  if (strstr(options,"CONTINUOUS"))
    Fsolver->dotMultiplyKernel(Fmesh->Nelements*Fmesh->Np, Fmesh->ogs->o_invDegree, o_x, o_x);

  precon->coarsenKernel(mesh->Nelements, o_R, o_x, o_Rx);

  if (strstr(options,"CONTINUOUS")) {
    //solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, o_Rx, o_Rx);
    ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, o_Rx, dfloatString, "add");  
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Rx);
  }
}

void prolongateQuad2D(void **args, occa::memory &o_x, occa::memory &o_Px) {

  solver_t *solver = (solver_t *) args[0];
  char *options    = (char *) args[1];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_R = solver->o_R;

  // if (strstr(options,"CONTINUOUS")) {
  //   if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_x);
  //   ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, o_x, dfloatString, "add");  
  //   solver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, o_x, o_x);
  // }

  precon->prolongateKernel(mesh->Nelements, o_R, o_x, o_Px);
}

void ellipticGather(void **args, occa::memory &o_x, occa::memory &o_Gx) {

  solver_t *solver = (solver_t *) args[0];
  ogs_t *ogs       = (ogs_t *) args[1];
  occa::memory *o_s= (occa::memory *) args[2];
  char *options    = (char *) args[3];
  mesh_t *mesh     = solver->mesh;

  meshParallelGather(mesh, ogs, o_x, o_Gx);  
  solver->dotMultiplyKernel(ogs->Ngather, ogs->o_gatherInvDegree, o_Gx, o_Gx);
}

void ellipticScatter(void **args, occa::memory &o_x, occa::memory &o_Sx) {

  solver_t *solver = (solver_t *) args[0];
  ogs_t *ogs       = (ogs_t *) args[1];
  occa::memory *o_s= (occa::memory *) args[2];
  char *options    = (char *) args[3];
  mesh_t *mesh     = solver->mesh;

  meshParallelScatter(mesh, ogs, o_x, o_Sx);  
}

void buildCoarsenerQuad2D(solver_t* solver, mesh2D **meshLevels, int Nf, int Nc, const char* options);

void ellipticMultiGridSetupQuad2D(solver_t *solver, precon_t* precon,
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
    meshLoadReferenceNodesQuad2D(meshLevels[n], n);
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
  } else if (strstr(options,"HALFDOFS")) {
    // pick the degrees so the dofs of each level halfs (roughly)
    //start by counting the number of levels neccessary
    numLevels = 1;
    int degree = mesh->N;
    int dofs = meshLevels[degree]->Np;
    while (dofs>4) {
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
    while (dofs>4) {
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

  //maually build multigrid levels
  precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);
  agmgLevel **levels = precon->parAlmond->levels;

  //build a solver struct for every degree
  solver_t **solversN = (solver_t**) calloc(mesh->N+1,sizeof(solver_t*));
  solversN[mesh->N] = solver; //top level
  for (int n=1;n<numLevels;n++) {  //build solver for this degree
    int Nf = levelDegree[n-1];
    int Nc = levelDegree[n];
    printf("=============BUIDLING MULTIGRID LEVEL OF DEGREE %d==================\n", Nc);
    solversN[Nc] = ellipticBuildMultigridLevelQuad2D(solver,Nc,Nf,BCType,options);
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
    levels[n]->device_Ax = AxQuad2D;

    levels[n]->smoothArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->smoothArgs[0] = (void *) solverL;
    levels[n]->smoothArgs[1] = (void *) levels[n];

    levels[n]->Nrows = mesh->Nelements*solverL->mesh->Np;
    levels[n]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*solverL->mesh->Np;

    if (strstr(options,"CHEBYSHEV")) {
      levels[n]->device_smooth = smoothChebyshevQuad2D;

      levels[n]->smootherResidual = (dfloat *) calloc(levels[n]->Ncols,sizeof(dfloat));

      // extra storage for smoothing op
      levels[n]->o_smootherResidual = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
      levels[n]->o_smootherResidual2 = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
      levels[n]->o_smootherUpdate = mesh->device.malloc(levels[n]->Ncols*sizeof(dfloat),levels[n]->smootherResidual);
    } else {
      levels[n]->device_smooth = smoothQuad2D;

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
  dlong nnzCoarseA;
  ogs_t *coarseogs;

  solver_t* solverL = solversN[1];

  hlong *coarseGlobalStarts = (hlong*) calloc(size+1, sizeof(hlong));

  if (strstr(options,"IPDG")) {
    ellipticBuildIpdgQuad2D(solverL->mesh, tau, lambda, BCType, &coarseA, &nnzCoarseA,coarseGlobalStarts, options);
  } else if (strstr(options,"CONTINUOUS")) {
    ellipticBuildContinuousQuad2D(solverL,lambda,&coarseA,&nnzCoarseA,&coarseogs,coarseGlobalStarts,options);
  }

  hlong *Rows = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  hlong *Cols = (hlong *) calloc(nnzCoarseA, sizeof(hlong));
  dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

  for (dlong i=0;i<nnzCoarseA;i++) {
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
    coarseLevel->o_Srhs = solverL->mesh->device.malloc(solverL->mesh->Np*solverL->mesh->Nelements*sizeof(dfloat),coarseLevel->Srhs);
    coarseLevel->o_Sx   = solverL->mesh->device.malloc(solverL->mesh->Np*solverL->mesh->Nelements*sizeof(dfloat),coarseLevel->Sx);

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
    solver_t *solverF = solversN[Nf];
    buildCoarsenerQuad2D(solverL, meshLevels, Nf, Nc, options);
    
    levels[n]->coarsenArgs = (void **) calloc(3,sizeof(void*));
    levels[n]->coarsenArgs[0] = (void *) solverL;
    levels[n]->coarsenArgs[1] = (void *) options;
    levels[n]->coarsenArgs[2] = (void *) solverF;

    levels[n]->prolongateArgs = levels[n]->coarsenArgs;
    
    levels[n]->device_coarsen = coarsenQuad2D;
    levels[n]->device_prolongate = prolongateQuad2D;
  }

  for (int n=1;n<mesh->N+1;n++) free(meshLevels[n]);
  free(meshLevels);
}

void buildCoarsenerQuad2D(solver_t* solver, mesh2D **meshLevels, int Nf, int Nc, const char* options) {

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
  solver->R = (dfloat *) calloc(NqFine*NqCoarse,sizeof(dfloat));
  for (int i=0;i<NqCoarse;i++) {
    for (int j=0;j<NqFine;j++) {
      solver->R[i*NqFine+j] = P[j*NqCoarse+i];
    }
  }
  solver->o_R = solver->mesh->device.malloc(NqFine*NqCoarse*sizeof(dfloat), solver->R);

  free(P); free(Ptmp);
}