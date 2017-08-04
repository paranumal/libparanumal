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
  occa::memory *o_V = (occa::memory *) args[1];

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  precon->coarsenKernel(mesh->Nelements, *o_V, o_x, o_Rx);
}

void prolongateTri2D(void **args, occa::memory &o_x, occa::memory &o_Px) {

  solver_t *solver = (solver_t *) args[0];
  occa::memory *o_V = (occa::memory *) args[1];

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  precon->prolongateKernel(mesh->Nelements, *o_V, o_x, o_Px);
}

dfloat *buildCoarsenerTri2D(mesh2D** meshLevels, int N, int Nc) {

  //TODO We can build the coarsen matrix either from the interRaise or interpLower matrices. Need to check which is better

  //use the Raise for now (essentally an L2 projection)

  int NpFine   = meshLevels[N]->Np;
  int NpCoarse = meshLevels[Nc]->Np;

  printf("allocating P as %d by %d\n", NpCoarse, NpFine);

  dfloat *P    = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
  dfloat *Ptmp = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));

  //initialize P as identity
  for (int i=0;i<NpCoarse;i++) P[i*NpCoarse+i] = 1.0;

  for (int n=Nc;n<N;n++) {

    int Npp1 = meshLevels[n+1]->Np;
    int Np   = meshLevels[n]->Np;

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

  //the coarsen matrix is P^T
  dfloat *R    = (dfloat *) calloc(NpFine*NpCoarse,sizeof(dfloat));
  for (int i=0;i<NpCoarse;i++) {
    for (int j=0;j<NpFine;j++) {
      R[i*NpFine+j] = P[j*NpCoarse+i];
    }
  }

  free(P); free(Ptmp);

  return R;
}

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon,
                                dfloat tau, dfloat lambda, iint *BCType,
                                const char *options, const char *parAlmondOptions) {

  void (*smoother)(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);

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

  int numLevels = mesh->N;
  int *levelDegree = (int *) calloc(numLevels,sizeof(int));
  for (int n=0;n<numLevels;n++) levelDegree[n] = mesh->N - n; //all degrees

  //maually build multigrid levels
  precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);

  //add top level
  precon->parAlmond->numLevels++;
  agmgLevel **levels = precon->parAlmond->levels;

  levels[0] = (agmgLevel *) calloc(1,sizeof(agmgLevel));

  dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
  *vlambda = lambda;
  levels[0]->AxArgs = (void **) calloc(3,sizeof(void*));
  levels[0]->AxArgs[0] = (void *) solver;
  levels[0]->AxArgs[1] = (void *) vlambda;
  levels[0]->AxArgs[2] = (void *) options;
  levels[0]->device_Ax = AxTri2D;

  levels[0]->smoothArgs = (void **) calloc(2,sizeof(void*));
  levels[0]->smoothArgs[0] = (void *) solver;
  levels[0]->smoothArgs[1] = (void *) levels[0];
  levels[0]->device_smooth = smoothTri2D;

  levels[0]->smootherArgs = (void **) calloc(1,sizeof(void*));
  levels[0]->smootherArgs[0] = (void *) solver;

  //set up the fine problem smoothing
  //TODO the weight should be found by estimating the max eigenvalue of Smoother*A, via Arnoldi?
  //  for now, the weights can be adjusted for stability here
  if (strstr(options, "IPDG")) {
    if(strstr(options, "OVERLAPPINGPATCH")){
      dfloat weight = 0.08; //stability weighting for smoother
      ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, tau, lambda, BCType, weight, options);
      levels[0]->device_smoother = overlappingPatchIpdg;

    } else if(strstr(options, "EXACTFULLPATCH")){
      dfloat weight = 0.85; //stability weighting for smoother
      ellipticSetupSmootherExactFullPatchIpdg(solver, precon, tau, lambda, BCType, weight, options);
      levels[0]->device_smoother = exactFullPatchIpdg;

    } else if(strstr(options, "APPROXFULLPATCH")){
      dfloat weight = 0.25; //stability weighting for smoother
      ellipticSetupSmootherApproxFullPatchIpdg(solver, precon, tau, lambda, BCType, weight, options);
      levels[0]->device_smoother = approxFullPatchIpdg;

    } else { //default to damped jacobi
      dfloat weight = 0.65; //stability weighting for smoother
      ellipticSetupSmootherDampedJacobiIpdg(solver, precon, tau, lambda, BCType, weight, options);
      levels[0]->device_smoother = dampedJacobi;
    }
  }

  levels[0]->Nrows = mesh->Nelements*mesh->Np;
  levels[0]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  // extra storage for smoothing op
  levels[0]->o_smootherResidual = mesh->device.malloc(levels[0]->Ncols*sizeof(dfloat),levels[0]->x);


  dfloat **V = (dfloat **) calloc(numLevels+1,sizeof(dfloat*));
  occa::memory *o_V = (occa::memory *) calloc(numLevels+1,sizeof(occa::memory));

  //add intermediate p levels
  for (int n=1;n<numLevels;n++) {
    //build ops for this level
    printf("=============BUIDLING MULTIGRID LEVEL OF DEGREE %d==================\n", levelDegree[n]);
    solver_t *solverL = ellipticBuildMultigridLevelTri2D(solver,levelDegree,n,options);

    //check if we're at the degree 1 problem
    if (n==numLevels-1) {
      // build degree 1 matrix problem
      nonZero_t *coarseA;
      iint nnzCoarseA;
      hgs_t *coarsehgs;
      dfloat *V1;

      iint *coarseGlobalStarts = (iint*) calloc(size+1, sizeof(iint));

      ellipticCoarsePreconditionerSetupTri2D(mesh, precon, tau, lambda, BCType,
                                             &V1, &coarseA, &nnzCoarseA,
                                             &coarsehgs, coarseGlobalStarts, options);

      precon->o_V1  = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), V1);

      iint *Rows = (iint *) calloc(nnzCoarseA, sizeof(iint));
      iint *Cols = (iint *) calloc(nnzCoarseA, sizeof(iint));
      dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

      for (iint n=0;n<nnzCoarseA;n++) {
        Rows[n] = coarseA[n].row;
        Cols[n] = coarseA[n].col;
        Vals[n] = coarseA[n].val;
      }

      // build amg starting at level 1
      parAlmondAgmgSetup(precon->parAlmond,
                         coarseGlobalStarts,
                         nnzCoarseA,
                         Rows,
                         Cols,
                         Vals,
                         coarsehgs);

      free(coarseA); free(Rows); free(Cols); free(Vals);
    } else {
      //build the level manually
      precon->parAlmond->numLevels++;
      levels[n] = (agmgLevel *) calloc(1,sizeof(agmgLevel));
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
    levels[n]->device_smooth = smoothTri2D;

    levels[n]->smootherArgs = (void **) calloc(1,sizeof(void*));
    levels[n]->smootherArgs[0] = (void *) solverL;

    //set up the fine problem smoothing
    //TODO the weight should be found by estimating the max eigenvalue of Smoother*A, via Arnoldi?
    //  for now, the weights can be adjusted for stability here
    if (strstr(options, "IPDG")) {
      if(strstr(options, "OVERLAPPINGPATCH")){
        dfloat weight = 0.08; //stability weighting for smoother
        ellipticSetupSmootherOverlappingPatchIpdg(solverL, solverL->precon, tau, lambda, BCType, weight, options);
        levels[n]->device_smoother = overlappingPatchIpdg;

      } else if(strstr(options, "EXACTFULLPATCH")){
        dfloat weight = 0.85; //stability weighting for smoother
        ellipticSetupSmootherExactFullPatchIpdg(solverL, solverL->precon, tau, lambda, BCType, weight, options);
        levels[n]->device_smoother = exactFullPatchIpdg;

      } else if(strstr(options, "APPROXFULLPATCH")){
        dfloat weight = 0.25; //stability weighting for smoother
        ellipticSetupSmootherApproxFullPatchIpdg(solverL, solverL->precon, tau, lambda, BCType, weight, options);
        levels[n]->device_smoother = approxFullPatchIpdg;

      } else { //default to damped jacobi
        dfloat weight = 0.65; //stability weighting for smoother
        ellipticSetupSmootherDampedJacobiIpdg(solverL, solverL->precon, tau, lambda, BCType, weight, options);
        levels[n]->device_smoother = dampedJacobi;
      }
    }

    levels[n]->Nrows = mesh->Nelements*solverL->mesh->Np;
    levels[n]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*solverL->mesh->Np;

    // extra storage for smoothing op
    levels[n]->o_smootherResidual = mesh->device.malloc(levels[0]->Ncols*sizeof(dfloat));

    //build coarsen and prologation ops
    int N = levelDegree[n-1];
    int Nc = levelDegree[n];

    printf("Nfine = %d, NCoarse = %d\n", N, Nc);

    V[n] = buildCoarsenerTri2D(meshLevels, N, Nc);

    printf("V dims = %d by %d\n",meshLevels[Nc]->Np, meshLevels[N]->Np);

    o_V[n] = mesh->device.malloc(meshLevels[N]->Np*meshLevels[Nc]->Np*sizeof(dfloat), V[n]);

    levels[n]->coarsenArgs = (void **) calloc(2,sizeof(void*));
    levels[n]->coarsenArgs[0] = (void *) solverL;
    levels[n]->coarsenArgs[1] = (void *) &(o_V[n]);
    levels[n]->device_coarsen = coarsenTri2D;

    levels[n]->prolongateArgs = levels[n]->coarsenArgs;
    levels[n]->device_prolongate = prolongateTri2D;
  }

  for (int n=1;n<mesh->N+1;n++) free(meshLevels[n]);
  free(meshLevels);
}