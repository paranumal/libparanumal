#include "elliptic.h"

void ellipticPreconditionerSetup(elliptic_t *elliptic, ogs_t *ogs, dfloat tau, dfloat lambda, int *BCType, const char *options, const char *parAlmondOptions){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = elliptic->mesh;
  precon_t *precon = elliptic->precon;

  if(strstr(options, "FULLALMOND")){ //build full A matrix and pass to Almond
    dlong nnz;
    nonZero_t *A;

    hlong *globalStarts = (hlong*) calloc(size+1, sizeof(hlong));

    int basisNp = mesh->Np;
    dfloat *basis = NULL;

    if (strstr(options,"BERN")) basis = mesh->VB;

    if (strstr(options,"IPDG")) {
      ellipticBuildIpdgTri2D(mesh, basisNp, basis, tau, lambda, BCType, &A, &nnz, globalStarts, options);
    } else if (strstr(options,"CONTINUOUS")) {
      ellipticBuildContinuousTri2D(elliptic,lambda,&A,&nnz, &(precon->ogs), globalStarts, options);
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
                       elliptic->allNeumann,
                       elliptic->allNeumannPenalty);
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
      baseLevel->gatherArgs[0] = (void *) elliptic;
      baseLevel->gatherArgs[1] = (void *) precon->ogs;
      baseLevel->gatherArgs[2] = (void *) &(baseLevel->o_Sx);
      baseLevel->gatherArgs[3] = (void *) options;
      baseLevel->scatterArgs = baseLevel->gatherArgs;

      baseLevel->device_gather  = ellipticGather;
      baseLevel->device_scatter = ellipticScatter;        
    }

    if (strstr(options,"MATRIXFREE")&&strstr(options,"IPDG")) { //swap the top AMG level ops for matrix free versions
      agmgLevel *baseLevel = precon->parAlmond->levels[0];

      dfloat *vlambda = (dfloat *) calloc(1,sizeof(dfloat));
      *vlambda = lambda;
      baseLevel->AxArgs = (void **) calloc(3,sizeof(void*));
      baseLevel->AxArgs[0] = (void *) elliptic;
      baseLevel->AxArgs[1] = (void *) vlambda;
      baseLevel->AxArgs[2] = (void *) options;
      baseLevel->device_Ax = AxTri2D;

      baseLevel->smoothArgs = (void **) calloc(2,sizeof(void*));
      baseLevel->smoothArgs[0] = (void *) elliptic;
      baseLevel->smoothArgs[1] = (void *) baseLevel;
      baseLevel->device_smooth = smoothTri2D;

      baseLevel->smootherArgs = (void **) calloc(1,sizeof(void*));
      baseLevel->smootherArgs[0] = (void *) elliptic;

      baseLevel->Nrows = mesh->Nelements*mesh->Np;
      baseLevel->Ncols = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

      // extra storage for smoothing op
      baseLevel->o_smootherResidual = mesh->device.malloc(baseLevel->Ncols*sizeof(dfloat),baseLevel->x);

      dfloat rateTolerance;    // 0 - accept not approximate patches, 1 - accept all approximate patches
      if(strstr(options, "EXACT")){
        rateTolerance = 0.0;
      } else {
        rateTolerance = 1.0;
      }

      //set up the fine problem smoothing
      if(strstr(options, "LOCALPATCH")){
        ellipticSetupSmootherLocalPatch(elliptic, precon, baseLevel, tau, lambda, BCType, rateTolerance, options);
      } else { //default to damped jacobi
        ellipticSetupSmootherDampedJacobi(elliptic, precon, baseLevel, tau, lambda, BCType, options);
      }
    }

  } else if (strstr(options, "MASSMATRIX")){

    // compute inverse mass matrix
    dfloat *dfMMinv = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
    double *MMinv = (double*) calloc(mesh->Np*mesh->Np, sizeof(double));
    int *ipiv = (int*) calloc(mesh->Np, sizeof(int));
    int lwork = mesh->Np*mesh->Np;
    double *work = (double*) calloc(lwork, sizeof(double));
    int info;
    for(int n=0;n<mesh->Np*mesh->Np;++n){
      MMinv[n] = mesh->MM[n];
    }

    dgetrf_ (&(mesh->Np), &(mesh->Np), MMinv, &(mesh->Np), ipiv, &info);
    dgetri_ (&(mesh->Np), MMinv, &(mesh->Np), ipiv, work, &lwork, &info);
    if(info)
      printf("dgetrf/dgetri reports info = %d when inverting the reference mass matrix\n", info);

    for(int n=0;n<mesh->Np*mesh->Np;++n){
      dfMMinv[n] = MMinv[n];
    }

    precon->o_invMM = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), dfMMinv);

    free(MMinv); free(ipiv); free(work); free(dfMMinv);

  } else if(strstr(options, "MULTIGRID")){

    ellipticMultiGridSetupTri2D(elliptic,precon,tau,lambda,BCType,options,parAlmondOptions);

  } else if(strstr(options, "SEMFEM")) {

    ellipticSEMFEMSetupTri2D(elliptic,precon,tau,lambda,BCType,options,parAlmondOptions);

  } else if(strstr(options,"JACOBI")) {

    dfloat *invDiagA;

    ellipticBuildJacobiTri2D(elliptic,mesh,mesh->Np,NULL,tau, lambda, BCType, &invDiagA,options);

    precon->o_invDiagA = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), invDiagA);

    free(invDiagA);
  }
}
