#include "ellipticQuad2D.h"


void ellipticPreconditionerSetupQuad2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, int *BCType, const char *options, const char *parAlmondOptions){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  if(strstr(options, "FULLALMOND")){ //build full A matrix and pass to Almond
    int nnz;
    nonZero_t *A;
    hgs_t *hgs;

    int Nnum = mesh->Np*mesh->Nelements;
    int *globalStarts = (int*) calloc(size+1, sizeof(int));

    if (strstr(options,"IPDG")) {
      ellipticBuildIpdgQuad2D(mesh, tau, lambda, BCType, &A, &nnz,globalStarts, options);
    } else if (strstr(options,"CONTINUOUS")) {
      ellipticBuildContinuousQuad2D(mesh,lambda,&A,&nnz,&hgs,globalStarts, options);
    }
    

    int *Rows = (int *) calloc(nnz, sizeof(int));
    int *Cols = (int *) calloc(nnz, sizeof(int));
    dfloat *Vals = (dfloat*) calloc(nnz,sizeof(dfloat));

    for (int n=0;n<nnz;n++) {
      Rows[n] = A[n].row;
      Cols[n] = A[n].col;
      Vals[n] = A[n].val;
    }

    precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);
    parAlmondAgmgSetup(precon->parAlmond,
                       globalStarts,
                       nnz,
                       Rows,
                       Cols,
                       Vals,
                       hgs);

    free(A); free(Rows); free(Cols); free(Vals);

    if (strstr(options,"MATRIXFREE")&&strstr(options,"IPDG")) { //swap the top AMG level ops for matrix free versions
      agmgLevel **levels = precon->parAlmond->levels;

      dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
      *vlambda = lambda;
      levels[0]->AxArgs = (void **) calloc(3,sizeof(void*));
      levels[0]->AxArgs[0] = (void *) solver;
      levels[0]->AxArgs[1] = (void *) vlambda;
      levels[0]->AxArgs[2] = (void *) options;
      levels[0]->device_Ax = AxQuad2D;

      levels[0]->smoothArgs = (void **) calloc(2,sizeof(void*));
      levels[0]->smoothArgs[0] = (void *) solver;
      levels[0]->smoothArgs[1] = (void *) levels[0];
      levels[0]->device_smooth = smoothQuad2D;

      levels[0]->smootherArgs = (void **) calloc(1,sizeof(void*));
      levels[0]->smootherArgs[0] = (void *) solver;

      levels[0]->Nrows = mesh->Nelements*mesh->Np;
      levels[0]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

      // extra storage for smoothing op
      levels[0]->o_smootherResidual = mesh->device.malloc(levels[0]->Ncols*sizeof(dfloat),levels[0]->x);

            dfloat rateTolerance;    // 0 - accept not approximate patches, 1 - accept all approximate patches
      if(strstr(options, "EXACT")){
        rateTolerance = 0.0;
      } else {
        rateTolerance = 1.0;
      }

      //set up the fine problem smoothing
      if(strstr(options, "OVERLAPPINGPATCH")){
        ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "FULLPATCH")){
        ellipticSetupSmootherFullPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, rateTolerance, options);
      } else if(strstr(options, "FACEPATCH")){
        ellipticSetupSmootherFacePatchIpdg(solver, precon, levels[0], tau, lambda, BCType, rateTolerance, options);
      } else if(strstr(options, "LOCALPATCH")){
        ellipticSetupSmootherLocalPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, rateTolerance, options);
      } else { //default to damped jacobi
        ellipticSetupSmootherDampedJacobiIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      }
    }

  } else if(strstr(options, "OAS")){

    //create a dummy parAlmond level to store the OAS smoother
    agmgLevel *OASLevel = (agmgLevel *) calloc(1,sizeof(agmgLevel));
    precon->OASLevel = OASLevel;
    precon->OASsmoothArgs = (void **) calloc(2,sizeof(void*));
    precon->OASsmoothArgs[0] = (void *) solver;
    precon->OASsmoothArgs[1] = (void *) precon->OASLevel;

    dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
    *vlambda = lambda;

    OASLevel->AxArgs = (void **) calloc(3,sizeof(void*));
    OASLevel->AxArgs[0] = (void *) solver;
    OASLevel->AxArgs[1] = (void *) vlambda;
    OASLevel->AxArgs[2] = (void *) options;
    OASLevel->device_Ax = AxQuad2D;

    OASLevel->smootherArgs = (void **) calloc(1,sizeof(void*));
    OASLevel->smootherArgs[0] = (void *) solver;

    OASLevel->Nrows = mesh->Nelements*mesh->Np;
    OASLevel->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

    // extra storage for smoothing op
    OASLevel->o_smootherResidual = mesh->device.malloc(OASLevel->Ncols*sizeof(dfloat));

    dfloat rateTolerance;    // 0 - accept not approximate patches, 1 - accept all approximate patches
    if(strstr(options, "EXACT")){
      rateTolerance = 0.0;
    } else {
      rateTolerance = 1.0;
    }

    //set up the fine problem smoothing
    if(strstr(options, "OVERLAPPINGPATCH")){
      ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, options);
    } else if(strstr(options, "FULLPATCH")){
      ellipticSetupSmootherFullPatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, rateTolerance, options);
    } else if(strstr(options, "FACEPATCH")){
      ellipticSetupSmootherFacePatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, rateTolerance, options);
    } else if(strstr(options, "LOCALPATCH")){
      ellipticSetupSmootherLocalPatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, rateTolerance, options);
    } else { //default to damped jacobi
      ellipticSetupSmootherDampedJacobiIpdg(solver, precon, OASLevel, tau, lambda, BCType, options);
    }


    // coarse grid preconditioner
    occaTimerTic(mesh->device,"CoarsePreconditionerSetup");
    nonZero_t *coarseA;
    int nnzCoarseA;
    hgs_t *coarsehgs;
    dfloat *V1;

    int *coarseGlobalStarts = (int*) calloc(size+1, sizeof(int));

    ellipticCoarsePreconditionerSetupQuad2D(mesh, precon, tau, lambda, BCType,
                                           &V1, &coarseA, &nnzCoarseA,
                                           &coarsehgs, coarseGlobalStarts, options);

    int Nnum = mesh->Nverts*(mesh->Nelements+mesh->totalHaloPairs);
    precon->o_V1  = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), V1);

    int *Rows = (int *) calloc(nnzCoarseA, sizeof(int));
    int *Cols = (int *) calloc(nnzCoarseA, sizeof(int));
    dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

    for (int n=0;n<nnzCoarseA;n++) {
      Rows[n] = coarseA[n].row;
      Cols[n] = coarseA[n].col;
      Vals[n] = coarseA[n].val;
    }

    precon->parAlmond = parAlmondInit(mesh, parAlmondOptions);
    parAlmondAgmgSetup(precon->parAlmond,
                       coarseGlobalStarts,
                       nnzCoarseA,
                       Rows,
                       Cols,
                       Vals,
                       coarsehgs);

    free(coarseA); free(Rows); free(Cols); free(Vals);

    precon->o_r1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->o_z1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->r1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    precon->z1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    occaTimerToc(mesh->device,"CoarsePreconditionerSetup");

  } else if(strstr(options, "MULTIGRID")){

    ellipticMultiGridSetupQuad2D(solver,precon,tau,lambda,BCType,options,parAlmondOptions);

  } else if(strstr(options,"JACOBI")) {

    dfloat *invDiagA;

    ellipticBuildJacobiIpdgQuad2D(mesh,mesh->Np,NULL,tau, lambda, BCType, &invDiagA,options);

    precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);
  }
}
