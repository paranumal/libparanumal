#include "ellipticTet3D.h"


void ellipticPreconditionerSetupTet3D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions){

  int rank, size;
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
      ellipticBuildIpdgTet3D(mesh, tau, lambda, BCType, &A, &nnz,globalStarts, options);
    } else if (strstr(options,"CONTINUOUS")) {
      ellipticBuildContinuousTet3D(mesh,lambda,&A,&nnz,&hgs,globalStarts, options);
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
      levels[0]->device_Ax = AxTet3D;

      levels[0]->smoothArgs = (void **) calloc(2,sizeof(void*));
      levels[0]->smoothArgs[0] = (void *) solver;
      levels[0]->smoothArgs[1] = (void *) levels[0];
      levels[0]->device_smooth = smoothTet3D;

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
        ellipticSetupSmootherApproxFullPatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "EXACTFACEPATCH")){
        ellipticSetupSmootherExactFacePatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "APPROXFACEPATCH")){
        ellipticSetupSmootherApproxFacePatchIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "EXACTBLOCKJACOBI")){
        ellipticSetupSmootherExactBlockJacobiIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else if(strstr(options, "APPROXBLOCKJACOBI")){
        ellipticSetupSmootherApproxBlockJacobiIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      } else { //default to damped jacobi
        ellipticSetupSmootherDampedJacobiIpdg(solver, precon, levels[0], tau, lambda, BCType, options);
      }
    }

  } else if (strstr(options, "MASSMATRIX")){

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

  } else if(strstr(options, "MULTIGRID")){

    ellipticMultiGridSetupTet3D(solver,precon,tau,lambda,BCType,options,parAlmondOptions);

  } else if(strstr(options,"JACOBI")) {

    dfloat *invDiagA;

    ellipticBuildJacobiIpdgTet3D(mesh,mesh->Np,NULL,tau, lambda, BCType, &invDiagA,options);

    precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);
  }
}