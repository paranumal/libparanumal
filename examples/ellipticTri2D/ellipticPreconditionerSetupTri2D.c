#include "ellipticTri2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

extern "C"
{
  void dgetrf_ (int *, int *, double *, int *, int *, int *);
  void dgetri_ (int *, double *, int *, int *, double *, int *, int *);
}

void ellipticBuildIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A, iint *nnzA,
                              hgs_t **hgs, iint *globalStarts, const char *options);

void ellipticBuildContinuousTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnz,
                              hgs_t **hgs, iint *globalStarts, const char* options);

void ellipticBuildPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, nonZero_t **A, iint *nnzA,
                                   hgs_t **hgs, iint *globalStarts,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA, dfloat **localA,
                                   const char *options);

void ellipticCoarsePreconditionerSetupTri2D(mesh_t *mesh, precon_t *precon, dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **V1, nonZero_t **A, iint *nnzA,
                                   hgs_t **hgs, iint *globalStarts, const char *options);

void ellipticBuildJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options);

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);


precon_t *ellipticPreconditionerSetupTri2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;
  precon_t *precon = (precon_t*) calloc(1, sizeof(precon_t));

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
    parAlmondAgmgSetup(precon->parAlmond,0,
                       globalStarts,
                       nnz,
                       Rows,
                       Cols,
                       Vals,
                       hgs,
                       parAlmondOptions);

    free(A); free(Rows); free(Cols); free(Vals);

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

    dfloat weight = 1.0; //stability weighting for smoother

    //set up the fine problem smoothing
    if (strstr(options, "IPDG")) {
      if(strstr(options, "OVERLAPPINGPATCH")){

        ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, tau, lambda, BCType, weight, options);

        precon->OASsmootherArgs = (void **) calloc(1,sizeof(void*));
        precon->OASsmootherArgs[0] = (void *) solver;
        precon->OASsmooth = overlappingPatchIpdg;

      } else if(strstr(options, "EXACTFULLPATCH")){

        ellipticSetupSmootherExactFullPatchIpdg(solver, precon, tau, lambda, BCType, weight, options);

        precon->OASsmootherArgs = (void **) calloc(1,sizeof(void*));
        precon->OASsmootherArgs[0] = (void *) solver;
        precon->OASsmooth = exactFullPatchIpdg;

      } else if(strstr(options, "APPROXFULLPATCH")){
        
        ellipticSetupSmootherApproxFullPatchIpdg(solver, precon, tau, lambda, BCType, weight, options);
        
        precon->OASsmootherArgs = (void **) calloc(1,sizeof(void*));
        precon->OASsmootherArgs[0] = (void *) solver;
        precon->OASsmooth = approxFullPatchIpdg;
      
      } else if(strstr(options, "DAMPEDJACOBI")){

        ellipticSetupSmootherDampedJacobiIpdg(solver, precon, tau, lambda, BCType, weight, options);
      
        precon->OASsmootherArgs = (void **) calloc(4,sizeof(void*));
        precon->OASsmootherArgs[0] = (void *) solver;
        precon->OASsmootherArgs[1] = (void *) &(precon->o_invDiagA);
        precon->OASsmootherArgs[2] = (void *) vlambda;
        precon->OASsmootherArgs[3] = (void *) options;

        precon->OASsmooth = dampedJacobi;
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
    parAlmondAgmgSetup(precon->parAlmond, 0,
                       coarseGlobalStarts,
                       nnzCoarseA,
                       Rows,
                       Cols,
                       Vals,
                       coarsehgs,
                       parAlmondOptions);

    free(coarseA); free(Rows); free(Cols); free(Vals);

    precon->o_r1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->o_z1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->r1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    precon->z1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    occaTimerToc(mesh->device,"CoarsePreconditionerSetup");

  } else if(strstr(options, "MULTIGRID")){

    ellipticMultiGridSetupTri2D(solver,precon,tau,lambda,BCType,options,parAlmondOptions);

  }


  return precon;
}
