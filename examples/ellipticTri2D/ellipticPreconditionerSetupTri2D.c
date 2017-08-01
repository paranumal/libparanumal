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

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon, dfloat tau, dfloat lambda, iint *BCType, const char *options);


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

    precon->parAlmond = parAlmondInit(mesh, options);
    parAlmondAgmgSetup(precon->parAlmond,0,
                       globalStarts,
                       nnz,
                       Rows,
                       Cols,
                       Vals,
                       hgs,
                       options);

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

    //set up the fine problem smoothing
    if (strstr(options, "IPDG")) {
      if(strstr(options, "OVERLAPPINGPATCH")){
        dfloat weight = 1.0; //stability weighting for smoother
        ellipticSetupSmootherOverlappingPatchIpdg(solver, precon, weight);
      }
    }


    if (strstr(options,"PATCHSOLVE")||strstr(options,"APPROXPATCH")
      ||strstr(options,"LOCALPATCH")) {
      iint nnz;
      nonZero_t *A;
      dfloat *invAP;
      dfloat *localInvA;
      iint Npatches;
      iint *patchesIndex;
      hgs_t *hgs;

      NpP = mesh->Np*(mesh->Nfaces+1);
      iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

      //initialize the full inverse operators on each 4 element patch
      ellipticBuildPatchesIpdgTri2D(mesh, mesh->Np, NULL, tau, lambda,
                                    BCType, &A, &nnz, &hgs, globalStarts,
                                    &Npatches, &patchesIndex, &invAP, &localInvA, options);

      Npatches = mesh->Nfaces*mesh->Nfaces*mesh->Nfaces+mesh->Nelements; // TW: HACK FOR THE MOMENT

      if (strstr(options,"PATCHSOLVE")) {
        precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);
        precon->o_patchesIndex = mesh->device.malloc(mesh->Nelements*sizeof(iint), patchesIndex);
      } else if (strstr(options,"APPROXPATCH")) {
        for (iint e =0;e<mesh->Nelements;e++) {
          //check if this is a boudary element
          int bcflag=0;
          for (int f=0;f<mesh->Nfaces;f++)
            if (mesh->EToE[e*mesh->Nfaces +f]<0) bcflag++;

          if (bcflag) break;

          iint offset = e*NpP*NpP;
          for (int n=0;n<NpP*NpP;n++)
            invAP[offset+n] = mesh->invAP[n];
        }

        //set reference patch inverse
        //precon->o_invAP = mesh->device.malloc(NpP*NpP*sizeof(dfloat), mesh->invAP);
        // TW: HACK FOR THE MOMENT (FIX LATER)
        precon->o_invAP = mesh->device.malloc(Npatches*NpP*NpP*sizeof(dfloat),invAP);

        //rotated node ids of neighbouring element
        mesh->o_rmapP = mesh->device.malloc(mesh->Np*mesh->Nfaces*sizeof(iint),mesh->rmapP);
        mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToF);
      } else if (strstr(options,"LOCALPATCH")) {
        //initialize the full inverse operators on each 4 element patch
        ellipticBuildPatchesIpdgTri2D(mesh, mesh->Np, NULL, tau, lambda,
                                      BCType, &A, &nnz, &hgs, globalStarts,
                                      &Npatches, &patchesIndex, &invAP, &localInvA, options);

        precon->o_invAP = mesh->device.malloc(mesh->Nelements*mesh->Np*mesh->Np*sizeof(dfloat),localInvA);
      }

      dfloat *invDegree = (dfloat*) calloc(mesh->Nelements,sizeof(dfloat));
      for (iint e=0;e<mesh->Nelements;e++) {
        for (int f=0;f<mesh->Nfaces;f++)
            invDegree[e] += (mesh->EToE[e*mesh->Nfaces +f]<0) ? 0 : 1; //overlap degree = # of neighbours
        invDegree[e] = 1./invDegree[e];
      }
      precon->o_invDegreeAP = mesh->device.malloc(mesh->Nelements*sizeof(dfloat),invDegree);
      free(invDegree);

      mesh->o_EToE = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToE);

      //set storage for larger patch
      free(solver->zP);
      solver->zP = (dfloat*) calloc(mesh->Nelements*NpP,  sizeof(dfloat));
      solver->o_zP.free();
      solver->o_zP = mesh->device.malloc(mesh->Nelements*NpP*sizeof(dfloat), solver->zP);

      //build a gather for the overlapping nodes

      // every degree of freedom has its own global id
      /* so find number of elements on each rank */
      iint *rankNelements = (iint*) calloc(size, sizeof(iint));
      iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
      MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
        rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
      //find offsets
      for(iint r=0;r<size;++r)
        rankStarts[r+1] = rankStarts[r]+rankNelements[r];

      iint *globalIds = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(iint));
      iint *globalOwners = (iint*) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np, sizeof(iint));

      //use the offsets to set the global ids
      for (iint e =0;e<mesh->Nelements;e++) {
        for (int n=0;n<mesh->Np;n++) {
          globalIds[e*mesh->Np + n] = n + (e + rankStarts[rank])*mesh->Np;
          globalOwners[e*mesh->Np +n] = rank;
        }
      }

      /* do a halo exchange of global node numbers */
      if (mesh->totalHaloPairs) {
        iint *idSendBuffer = (iint *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(iint));
        meshHaloExchange(mesh, mesh->Np*sizeof(iint), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
        meshHaloExchange(mesh, mesh->Np*sizeof(iint), globalOwners, idSendBuffer, globalOwners + mesh->Nelements*mesh->Np);
        free(idSendBuffer);
      }

      iint *globalIdsP = (iint *) calloc(mesh->Nelements*NpP,sizeof(iint));
      iint *globalOwnersP = (iint*) calloc(mesh->Nelements*NpP, sizeof(iint));

      //use the global ids of the nodes to set the ids of each patch
      for (iint e =0;e<mesh->Nelements;e++) {
        for (int f=0;f<mesh->Nfaces;f++) {
          for (int n=0;n<mesh->Np;n++) {
            int eP = (f==0) ? e: mesh->EToE[e*mesh->Nfaces+f-1]; // load element e first
            int fP = (f==0) ? 0: mesh->EToF[e*mesh->Nfaces+f-1];

            if (eP>-1) {
              if (strstr(options,"APPROXPATCH")) {
                int id = mesh->rmapP[fP*mesh->Np+n];
                globalIdsP[e*NpP + f*mesh->Np + n] = globalIds[eP*mesh->Np + id];
                globalOwnersP[e*NpP + f*mesh->Np + n] = globalOwners[eP*mesh->Np + id];
              } else {
                globalIdsP[e*NpP + f*mesh->Np + n] = globalIds[eP*mesh->Np + n];
                globalOwnersP[e*NpP + f*mesh->Np + n] = globalOwners[eP*mesh->Np + n];
              }
            } else {
              globalIdsP[e*NpP + f*mesh->Np + n] = -1;
              globalOwnersP[e*NpP + f*mesh->Np + n] = -1;
            }
          }
        }
      }
      free(globalIds); free(globalOwners);
      free(rankNelements); free(rankStarts);

      //use the ordering to define a gather+scatter for assembly
      precon->hgsDg = meshParallelGatherSetup(mesh, NpP*mesh->Nelements, globalIdsP, globalOwnersP);
    }

    // coarse grid preconditioner (only continous elements)
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

    parAlmondAgmgSetup(precon->parAlmond, 0,
                       coarseGlobalStarts,
                       nnzCoarseA,
                       Rows,
                       Cols,
                       Vals,
                       coarsehgs,
                       options);

    free(coarseA); free(Rows); free(Cols); free(Vals);

    precon->o_r1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->o_z1 = mesh->device.malloc(Nnum*sizeof(dfloat));
    precon->r1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    precon->z1 = (dfloat*) malloc(Nnum*sizeof(dfloat));
    occaTimerToc(mesh->device,"CoarsePreconditionerSetup");

  } else if(strstr(options, "OMS")){

    ellipticMultiGridSetupTri2D(solver,precon,tau,lambda,BCType,options);

  }


  return precon;
}
