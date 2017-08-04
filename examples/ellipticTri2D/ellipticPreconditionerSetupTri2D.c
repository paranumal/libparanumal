#include "ellipticTri2D.h"


void ellipticPreconditionerSetupTri2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  if(strstr(options, "FULLALMOND")){ //build full A matrix and pass to Almond
    iint nnz;
    nonZero_t *A;
    hgs_t *hgs;

    iint Nnum = mesh->Np*mesh->Nelements;
    iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

    if (strstr(options,"IPDG")) {
      ellipticBuildIpdgTri2D(mesh, tau, lambda, BCType, &A, &nnz,&hgs,globalStarts, options);
    } else if (strstr(options,"CONTINUOUS")) {
      ellipticBuildContinuousTri2D(mesh,lambda,&A,&nnz,&hgs,globalStarts, options);
    }

    iint *Rows = (iint *) calloc(nnz, sizeof(iint));
    iint *Cols = (iint *) calloc(nnz, sizeof(iint));
    dfloat *Vals = (dfloat*) calloc(nnz,sizeof(dfloat));

    for (iint n=0;n<nnz;n++) {
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
      levels[0]->device_Ax = AxTri2D;

      levels[0]->smoothArgs = (void **) calloc(2,sizeof(void*));
      levels[0]->smoothArgs[0] = (void *) solver;
      levels[0]->smoothArgs[1] = (void *) levels[0];
      levels[0]->device_smooth = smoothTri2D;

      levels[0]->smootherArgs = (void **) calloc(1,sizeof(void*));
      levels[0]->smootherArgs[0] = (void *) solver;

      levels[0]->Nrows = mesh->Nelements*mesh->Np;
      levels[0]->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

      // extra storage for smoothing op
      levels[0]->o_smootherResidual = mesh->device.malloc(levels[0]->Ncols*sizeof(dfloat),levels[0]->x);

      //set up the fine problem smoothing
      if(strstr(options, "OVERLAPPINGPATCH")){
        ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "EXACTFULLPATCH")){
        ellipticSetupSmootherExactFullPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "APPROXFULLPATCH")){
        ellipticSetupSmootherApproxFullPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);\
      } else { //default to damped jacobi
        ellipticSetupSmootherDampedJacobiIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      }
    }

  } else if (strstr(options, "BLOCKJACOBI")){

    // compute inverse mass matrix
    dfloat *dfMMinv = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
    double *MMinv = (double*) calloc(mesh->Np*mesh->Np, sizeof(double));
    iint *ipiv = (iint*) calloc(mesh->Np, sizeof(iint));
    int lwork = mesh->Np*mesh->Np;
    double *work = (double*) calloc(lwork, sizeof(double));
    iint info;
    for(iint n=0;n<mesh->Np*mesh->Np;++n){
      MMinv[n] = mesh->MM[n];
    }

    dgetrf_ (&(mesh->Np), &(mesh->Np), MMinv, &(mesh->Np), ipiv, &info);
    dgetri_ (&(mesh->Np), MMinv, &(mesh->Np), ipiv, work, &lwork, &info);
    if(info)
      printf("dgetrf/dgetri reports info = %d when inverting the reference mass matrix\n", info);

    for(iint n=0;n<mesh->Np*mesh->Np;++n){
      dfMMinv[n] = MMinv[n];
    }

    precon->o_invMM = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), dfMMinv);

    free(MMinv); free(ipiv); free(work); free(dfMMinv);

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
    OASLevel->device_Ax = AxTri2D;

    OASLevel->smootherArgs = (void **) calloc(1,sizeof(void*));
    OASLevel->smootherArgs[0] = (void *) solver;

    OASLevel->Nrows = mesh->Nelements*mesh->Np;
    OASLevel->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

    // extra storage for smoothing op
    OASLevel->o_smootherResidual = mesh->device.malloc(OASLevel->Ncols*sizeof(dfloat));

    //set up the fine problem smoothing
    if (strstr(options, "IPDG")) {
      if(strstr(options, "OVERLAPPINGPATCH")){
        ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, options);
      } else if(strstr(options, "EXACTFULLPATCH")){
        ellipticSetupSmootherExactFullPatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, options);
      } else if(strstr(options, "APPROXFULLPATCH")){
        ellipticSetupSmootherApproxFullPatchIpdg(solver, precon, OASLevel, tau, lambda, BCType, options);
      } else { //default to damped jacobi
        ellipticSetupSmootherDampedJacobiIpdg(solver, precon, OASLevel, tau, lambda, BCType, options);
      }
    }


    // coarse grid preconditioner
    occaTimerTic(mesh->device,"CoarsePreconditionerSetup");
    nonZero_t *coarseA;
    iint nnzCoarseA;
    hgs_t *coarsehgs;
    dfloat *V1;

    iint *coarseGlobalStarts = (iint*) calloc(size+1, sizeof(iint));

    ellipticCoarsePreconditionerSetupTri2D(mesh, precon, tau, lambda, BCType,
                                           &V1, &coarseA, &nnzCoarseA,
                                           &coarsehgs, coarseGlobalStarts, options);

    iint Nnum = mesh->Nverts*(mesh->Nelements+mesh->totalHaloPairs);
    precon->o_V1  = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), V1);

    iint *Rows = (iint *) calloc(nnzCoarseA, sizeof(iint));
    iint *Cols = (iint *) calloc(nnzCoarseA, sizeof(iint));
    dfloat *Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

    for (iint n=0;n<nnzCoarseA;n++) {
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

    ellipticMultiGridSetupTri2D(solver,precon,tau,lambda,BCType,options,parAlmondOptions);

  }
}
