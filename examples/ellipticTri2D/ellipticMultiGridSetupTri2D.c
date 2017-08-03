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

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon,
                                dfloat tau, dfloat lambda, iint *BCType,
                                const char *options, const char *parAlmondOptions) {

  void (*smoother)(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;

  //maually build multigrid levels
  precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);

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

    } else if(strstr(options, "DAMPEDJACOBI")){
      dfloat weight = 0.65; //stability weighting for smoother
      ellipticSetupSmootherDampedJacobiIpdg(solver, precon, tau, lambda, BCType, weight, options);
      levels[0]->device_smoother = dampedJacobi;
    }
  }

  levels[0]->Nrows = mesh->Nelements*mesh->Np;
  levels[0]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  // extra storage for smoothing op
  levels[0]->o_smootherResidual = mesh->device.malloc(levels[0]->Ncols*sizeof(dfloat),levels[0]->x);

  // coarse grid preconditioner
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
  parAlmondAgmgSetup(precon->parAlmond, 1,
                     coarseGlobalStarts,
                     nnzCoarseA,
                     Rows,
                     Cols,
                     Vals,
                     coarsehgs);

  levels[1]->coarsenArgs = (void **) calloc(2,sizeof(void*));
  levels[1]->coarsenArgs[0] = (void *) solver;
  levels[1]->coarsenArgs[1] = (void *) &(precon->o_V1);
  levels[1]->device_coarsen = coarsenTri2D;

  levels[1]->prolongateArgs = levels[1]->coarsenArgs;
  levels[1]->device_prolongate = prolongateTri2D;


  free(coarseA); free(Rows); free(Cols); free(Vals);
}