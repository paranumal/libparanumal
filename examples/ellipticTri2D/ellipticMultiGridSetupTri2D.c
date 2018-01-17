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

  solver_t **solversN = (solver_t **) args[0];
  int *Nstart = (int *) args[1];
  int *Nend = (int *) args[2];
  char *options    = (char *) args[3];

  for (int n=*Nstart;n>=*Nend;n--) {
    solver_t *solver = solversN[n];
    mesh_t *mesh = solver->mesh;
    precon_t *precon = solver->precon;
    occa::memory o_R = solver->o_R;

    //printf("coarsening degree %d to degree %d \n", n+1, n);

    occa::memory o_y, o_Ry;
    if (n==*Nstart) {
      o_y = o_x;
    } else {
      o_y = solversN[n+1]->o_Ry;
    }
    if (n==*Nend) {
      o_Ry = o_Rx;
    } else {
      o_Ry = solversN[n]->o_Ry;
    }

    precon->coarsenKernel(mesh->Nelements, o_R, o_y, o_Ry);

    // gather-scatter
    //sign correction for gs
    if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Ry, mesh->o_mapSgn, o_Ry);
    ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_Ry, o_Ry, dfloatString, "add");  
    if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Ry, mesh->o_mapSgn, o_Ry);       
    
    solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Ry, solver->o_invDegree, o_Ry);
  }
}

void prolongateTri2D(void **args, occa::memory &o_x, occa::memory &o_Px) {

  solver_t **solversN = (solver_t **) args[0];
  int *Nstart = (int *) args[1];
  int *Nend = (int *) args[2];
  char *options    = (char *) args[3];

  for (int n=*Nend;n<=*Nstart;n++) {
    solver_t *solver = solversN[n];
    mesh_t *mesh = solver->mesh;
    precon_t *precon = solver->precon;
    occa::memory o_R = solver->o_R;

    //printf("prologating degree %d to degree %d \n", n, n+1);

    occa::memory o_y, o_Py;
    if (n==*Nend) {
      o_y = o_x;
    } else {
      o_y = solversN[n]->o_Ry;
    }
    o_Py = solversN[n+1]->o_Ry;

    precon->prolongateKernel(mesh->Nelements, o_R, o_y, o_Py);

    // gather-scatter
    // solver = solversN[n+1];
    // mesh = solver->mesh;
    // //sign correction for gs
    // if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Py, mesh->o_mapSgn, o_Py);
    // ellipticParallelGatherScatterTri2D(mesh, solver->ogs, o_Py, o_Py, dfloatString, "add");  
    // if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Py, mesh->o_mapSgn, o_Py);       
    
    // solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Py, solver->o_invDegree, o_Py);
  }

  //add into Px
  ellipticScaledAdd(solversN[*Nstart+1], 1.0, solversN[*Nstart+1]->o_Ry, 1.0, o_Px);
}

void ellipticGather(void **args, occa::memory &o_x, occa::memory &o_Gx) {

  solver_t *solver = (solver_t *) args[0];
  hgs_t *hgs       = (hgs_t *) args[1];
  char *options    = (char *) args[2];
  mesh_t *mesh     = solver->mesh;

  if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_x, mesh->o_mapSgn, o_x);
  meshParallelGather(mesh, hgs, o_x, o_Gx);
  solver->dotMultiplyKernel(hgs->Ngather, hgs->o_invDegree, o_Gx, o_Gx);
}

void ellipticScatter(void **args, occa::memory &o_x, occa::memory &o_Sx) {

  solver_t *solver = (solver_t *) args[0];
  hgs_t *hgs       = (hgs_t *) args[1];
  char *options    = (char *) args[2];
  mesh_t *mesh     = solver->mesh;

  meshParallelScatter(mesh, hgs, o_x, o_Sx);
  if (strstr(options,"SPARSE")) solver->dotMultiplyKernel(mesh->Np*mesh->Nelements, o_Sx, mesh->o_mapSgn, o_Sx);
}

int buildCoarsenerTri2D(solver_t** solversN, int Nf, int Nc, const char* options);

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon,
                                dfloat tau, dfloat lambda, iint *BCType,
                                const char *options, const char *parAlmondOptions) {

  iint rank, size;
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

  //initialize parAlmond
  precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);
  agmgLevel **levels = precon->parAlmond->levels;

  //build a solver struct for every degree
  solver_t **solversN = (solver_t**) calloc(mesh->N+1,sizeof(solver_t*));
  solversN[mesh->N] = solver; //top level
  for (int n=mesh->N-1;n>0;n--) {  //build solver for this degree
    printf("=============BUIDLING MULTIGRID LEVEL OF DEGREE %d==================\n", n);
    solversN[n] = ellipticBuildMultigridLevelTri2D(solver,n,BCType,options);
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

    levels[n]->smootherArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->smootherArgs[0] = (void *) solverL;
    levels[n]->smootherArgs[1] = (void *) options;

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

      iint Nrows = levels[lev]->Nrows;

      iint minNrows=0, maxNrows=0, totalNrows=0;
      dfloat avgNrows;
      MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_IINT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&Nrows, &totalNrows, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
      avgNrows = (dfloat) totalNrows/size;

      if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
      MPI_Allreduce(&Nrows, &minNrows, 1, MPI_IINT, MPI_MIN, MPI_COMM_WORLD);

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
  iint nnzCoarseA;
  hgs_t *coarsehgs;

  solver_t* solverL = solversN[1];
  int basisNp = solverL->mesh->Np;
  dfloat *basis = NULL;

  if (strstr(options,"BERN")) basis = solverL->mesh->VB;

  iint *coarseGlobalStarts = (iint*) calloc(size+1, sizeof(iint));

  if (strstr(options,"IPDG")) {
    ellipticBuildIpdgTri2D(solverL->mesh, basisNp, basis, tau, lambda, BCType, &coarseA, &nnzCoarseA,coarseGlobalStarts, options);
  } else if (strstr(options,"BRDG")) {
    ellipticBuildBRdgTri2D(solverL->mesh, basisNp, basis, tau, lambda, BCType, &coarseA, &nnzCoarseA,coarseGlobalStarts, options);
  } else if (strstr(options,"CONTINUOUS")) {
    ellipticBuildContinuousTri2D(solverL->mesh,lambda,&coarseA,&nnzCoarseA,&coarsehgs,coarseGlobalStarts, options);
  }

  iint *Rows = (iint *) calloc(nnzCoarseA, sizeof(iint));
  iint *Cols = (iint *) calloc(nnzCoarseA, sizeof(iint));
  dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

  for (iint i=0;i<nnzCoarseA;i++) {
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
    coarseLevel->o_Srhs = mesh->device.malloc(solverL->mesh->Np*solverL->mesh->Nelements*sizeof(dfloat),coarseLevel->Srhs);
    coarseLevel->o_Sx   = mesh->device.malloc(solverL->mesh->Np*solverL->mesh->Nelements*sizeof(dfloat),coarseLevel->Sx);

    coarseLevel->gatherArgs = (void **) calloc(3,sizeof(void*));  
    coarseLevel->gatherArgs[0] = (void *) solverL;
    coarseLevel->gatherArgs[1] = (void *) coarsehgs;
    coarseLevel->gatherArgs[2] = (void *) options;
    coarseLevel->scatterArgs = coarseLevel->gatherArgs;

    coarseLevel->device_gather  = ellipticGather;
    coarseLevel->device_scatter = ellipticScatter;        
  }

  /* build coarsening and prologation operators to connect levels */
  int *NcoarseStart = (int *) calloc(numLevels+1,sizeof(int));
  int *NcoarseEnd = (int *) calloc(numLevels+1,sizeof(int));
  for(int n=1; n<numLevels; n++) {
    //build coarsen and prologation ops
    int N = levelDegree[n-1];
    int Nc = levelDegree[n];

    int NcoarsenOps = buildCoarsenerTri2D(solversN, N, Nc, options);
    NcoarseEnd[n] = Nc;
    NcoarseStart[n] = Nc+NcoarsenOps-1;
    
    levels[n]->coarsenArgs = (void **) calloc(4,sizeof(void*));
    levels[n]->coarsenArgs[0] = (void *) solversN;
    levels[n]->coarsenArgs[1] = (void *) (NcoarseStart+n);
    levels[n]->coarsenArgs[2] = (void *) (NcoarseEnd+n);
    levels[n]->coarsenArgs[3] = (void *) options;

    levels[n]->prolongateArgs = levels[n]->coarsenArgs;
    
    levels[n]->device_coarsen = coarsenTri2D;
    levels[n]->device_prolongate = prolongateTri2D;
  }

  for (int n=1;n<mesh->N+1;n++) free(meshLevels[n]);
  free(meshLevels);
}



int buildCoarsenerTri2D(solver_t** solversN, int Nf, int Nc, const char* options) {

  //use the Raise for now (essentally an L2 projection)

  if (strstr(options,"IPDG")) {
    int NpFine   = solversN[Nf]->mesh->Np;
    int NpCoarse = solversN[Nc]->mesh->Np;
    dfloat *P    = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
    dfloat *Ptmp = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));

    //initialize P as identity
    for (int i=0;i<NpCoarse;i++) P[i*NpCoarse+i] = 1.0;

    for (int n=Nc;n<Nf;n++) {

      int Npp1 = solversN[n+1]->mesh->Np;
      int Np   = solversN[n  ]->mesh->Np;

      //copy P
      for (int i=0;i<Np*NpCoarse;i++) Ptmp[i] = P[i];

      //Multiply by the raise op
      for (int i=0;i<Npp1;i++) {
        for (int j=0;j<NpCoarse;j++) {
          P[i*NpCoarse + j] = 0.;
          for (int k=0;k<Np;k++) {
            P[i*NpCoarse + j] += solversN[n]->mesh->interpRaise[i*Np+k]*Ptmp[k*NpCoarse + j];
          }
        }
      }
    }

    if (strstr(options,"BERN")) {
      dfloat* BBP = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
      for (iint j=0;j<NpFine;j++) {
        for (iint i=0;i<NpCoarse;i++) {
          for (iint k=0;k<NpCoarse;k++) {
            for (iint l=0;l<NpFine;l++) {
              BBP[i+j*NpCoarse] += solversN[Nf]->mesh->invVB[l+j*NpFine]*P[k+l*NpCoarse]*solversN[Nc]->mesh->VB[i+k*NpCoarse];
            }
          }
        }
      }
      for (iint j=0;j<NpFine;j++) {
        for (iint i=0;i<NpCoarse;i++) {
          P[i+j*NpCoarse] = BBP[i+j*NpCoarse];
        }
      }
      free(BBP);
    }

    //the coarsen matrix is P^T
    solversN[Nc]->R = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
    for (int i=0;i<NpCoarse;i++) {
      for (int j=0;j<NpFine;j++) {
        solversN[Nc]->R[i*NpFine+j] = P[j*NpCoarse+i];
      }
    }
    solversN[Nc]->o_R = solversN[Nc]->mesh->device.malloc(NpFine*NpCoarse*sizeof(dfloat), solversN[Nc]->R);

    free(P); free(Ptmp);

    //sizes for the coarsen and prolongation kernels. degree NFine to degree N
    occa::kernelInfo kernelInfo = solversN[Nc]->kernelInfo;

    kernelInfo.addDefine("p_NpFine", NpFine);
    kernelInfo.addDefine("p_NpCoarse", NpCoarse);

    solversN[Nc]->precon->coarsenKernel =
      solversN[Nc]->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
               "ellipticPreconCoarsen",
               kernelInfo);

    solversN[Nc]->precon->prolongateKernel =
      solversN[Nc]->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
               "ellipticPreconProlongate",
               kernelInfo);
    
    return 1; //only 1 coarsening required

  } else {
    int NcoarsenOps = Nf-Nc;

    solversN[Nf]->Ry = (dfloat*) calloc(solversN[Nf]->mesh->Nelements*solversN[Nf]->mesh->Np,sizeof(dfloat));
    solversN[Nf]->o_Ry = solversN[Nf]->mesh->device.malloc(solversN[Nf]->mesh->Nelements*solversN[Nf]->mesh->Np*sizeof(dfloat),solversN[Nf]->Ry);
    for (int n=Nc;n<Nf;n++) {
      int NpFine   = solversN[n+1]->mesh->Np;
      int NpCoarse = solversN[n  ]->mesh->Np;

      solversN[n]->R = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
      for (int i=0;i<NpCoarse;i++) {
        for (int j=0;j<NpFine;j++) {
          solversN[n]->R[i*NpFine+j] = solversN[n]->mesh->interpRaise[j*NpCoarse+i];
        }
      }
      solversN[n]->o_R  = solversN[n]->mesh->device.malloc(NpFine*NpCoarse*sizeof(dfloat), solversN[n]->R);
      solversN[n]->Ry   = (dfloat*) calloc(solversN[n]->mesh->Nelements*solversN[n]->mesh->Np,sizeof(dfloat));
      solversN[n]->o_Ry = solversN[n]->mesh->device.malloc(solversN[n]->mesh->Nelements*solversN[n]->mesh->Np*sizeof(dfloat),solversN[n]->Ry);
    
      //sizes for the coarsen and prolongation kernels. degree NFine to degree N
      occa::kernelInfo kernelInfo = solversN[n]->kernelInfo;
      
      kernelInfo.addDefine("p_NpFine", NpFine);
      kernelInfo.addDefine("p_NpCoarse", NpCoarse);

      solversN[n]->precon->coarsenKernel =
        solversN[n]->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconCoarsen.okl",
                 "ellipticPreconCoarsen",
                 kernelInfo);

      solversN[n]->precon->prolongateKernel =
        solversN[n]->mesh->device.buildKernelFromSource(DHOLMES "/okl/ellipticPreconProlongate.okl",
                 "ellipticPreconProlongate",
                 kernelInfo);
    }

    return NcoarsenOps;
  }
}