#include "ellipticHex3D.h"


void ellipticPreconditionerSetupHex3D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, int *BCType, const char *options, const char *parAlmondOptions){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh3D *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  if(strstr(options, "FULLALMOND")){ //build full A matrix and pass to Almond
    dlong nnz;
    nonZero_t *A;
    
    hlong *globalStarts = (hlong*) calloc(size+1, sizeof(hlong));

    if (strstr(options,"IPDG")) {
      ellipticBuildIpdgHex3D(mesh, tau, lambda, BCType, &A, &nnz,globalStarts, options);
    } else if (strstr(options,"CONTINUOUS")) {
      ellipticBuildContinuousHex3D(solver,lambda,&A,&nnz,&(precon->ogs),globalStarts, options);
    }
    

    hlong *Rows = (hlong *) calloc(nnz, sizeof(hlong));
    hlong *Cols = (hlong *) calloc(nnz, sizeof(hlong));
    dfloat *Vals = (dfloat*) calloc(nnz,sizeof(dfloat));

    for (dlong n=0;n<nnz;n++) {
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
                       solver->allNeumann,
                       solver->allNeumannPenalty);
    free(A); free(Rows); free(Cols); free(Vals);

    if (strstr(options,"CONTINUOUS")) {//tell parAlmond to gather this level
      agmgLevel *baseLevel = precon->parAlmond->levels[0];

      baseLevel->gatherLevel = true;
      baseLevel->Srhs = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
      baseLevel->Sx   = (dfloat*) calloc(mesh->Np*mesh->Nelements,sizeof(dfloat));
      baseLevel->o_Srhs = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat));
      baseLevel->o_Sx   = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat));

      baseLevel->weightedInnerProds = false;

      baseLevel->gatherArgs = (void **) calloc(4,sizeof(void*));  
      baseLevel->gatherArgs[0] = (void *) solver;
      baseLevel->gatherArgs[1] = (void *) precon->ogs;
      baseLevel->gatherArgs[2] = (void *) &(baseLevel->o_Sx);
      baseLevel->gatherArgs[3] = (void *) options;
      baseLevel->scatterArgs = baseLevel->gatherArgs;

      baseLevel->device_gather  = ellipticGather;
      baseLevel->device_scatter = ellipticScatter;        
    }

    if (strstr(options,"MATRIXFREE")&&strstr(options,"IPDG")) { //swap the top AMG level ops for matrix free versions
      agmgLevel **levels = precon->parAlmond->levels;

      dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
      *vlambda = lambda;
      levels[0]->AxArgs = (void **) calloc(3,sizeof(void*));
      levels[0]->AxArgs[0] = (void *) solver;
      levels[0]->AxArgs[1] = (void *) vlambda;
      levels[0]->AxArgs[2] = (void *) options;
      levels[0]->device_Ax = AxHex3D;

      levels[0]->smoothArgs = (void **) calloc(2,sizeof(void*));
      levels[0]->smoothArgs[0] = (void *) solver;
      levels[0]->smoothArgs[1] = (void *) levels[0];
      levels[0]->device_smooth = smoothHex3D;

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
        ellipticSetupSmootherOverlappingPatch(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "LOCALPATCH")){
        ellipticSetupSmootherLocalPatch(solver, precon, levels[0], tau, lambda, BCType, rateTolerance, options);
      } else { //default to damped jacobi
        ellipticSetupSmootherDampedJacobi(solver, precon, levels[0], tau, lambda, BCType, options);
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
    OASLevel->device_Ax = AxHex3D;

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
      ellipticSetupSmootherOverlappingPatch(solver, precon, OASLevel, tau, lambda, BCType, options);
    } else if(strstr(options, "LOCALPATCH")){
      ellipticSetupSmootherLocalPatch(solver, precon, OASLevel, tau, lambda, BCType, rateTolerance, options);
    } else { //default to damped jacobi
      ellipticSetupSmootherDampedJacobi(solver, precon, OASLevel, tau, lambda, BCType, options);
    }


    // coarse grid preconditioner
    occaTimerTic(mesh->device,"CoarsePreconditionerSetup");
    //TODO need to fix this to use the normal builder for OAS to work again
    /*
    nonZero_t *coarseA;
    int nnzCoarseA;
    hgs_t *coarsehgs;
    dfloat *V1;

    int *coarseGlobalStarts = (int*) calloc(size+1, sizeof(int));

    ellipticCoarsePreconditionerSetupHex3D(mesh, precon, tau, lambda, BCType,
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
    */
  } else if(strstr(options, "MULTIGRID")){

    ellipticMultiGridSetupHex3D(solver,precon,tau,lambda,BCType,options,parAlmondOptions);

  } else if(strstr(options, "SEMFEM")) {

    ellipticSEMFEMSetupHex3D(solver,precon,tau,lambda,BCType,options,parAlmondOptions);

  } else if(strstr(options,"JACOBI")) {

    dfloat *invDiagA;

    ellipticBuildJacobiHex3D(solver,tau, lambda, BCType, &invDiagA,options);

    precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);

    free(invDiagA);    
  }
}
