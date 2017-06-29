#include "ellipticTri2D.h"

typedef struct{

  iint localId;
  iint baseId;
  iint haloFlag;

} preconGatherInfo_t;

int parallelCompareBaseId(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;

}

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

// compare on global indices
int parallelCompareRowColumn(const void *a, const void *b);

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
                                   dfloat **patchesInvA, dfloat **localA,
                                   const char *options);

precon_t *ellipticPreconditionerSetupTri2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options){

  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh2D *mesh = solver->mesh;

  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;

  precon_t *precon = (precon_t*) calloc(1, sizeof(precon_t));

  if(strstr(options, "OAS")||strstr(options, "OMS")){
    //set up the fine problem smoothing

    // build gather-scatter
    iint NpP = mesh->Np + mesh->Nfaces*mesh->Nfp;

    // build gather-scatter for overlapping patches
    iint *allNelements = (iint*) calloc(size, sizeof(iint));
    MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
                allNelements, 1, MPI_IINT, MPI_COMM_WORLD);

    // offsets
    iint *startElement = (iint*) calloc(size, sizeof(iint));
    for(iint r=1;r<size;++r){
      startElement[r] = startElement[r-1]+allNelements[r-1];
    }

    // 1-indexed numbering of nodes on this process
    iint *localNums = (iint*) calloc((Nlocal+Nhalo), sizeof(iint));
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint n=0;n<mesh->Np;++n){
        localNums[e*mesh->Np+n] = 1 + e*mesh->Np + n + startElement[rank]*mesh->Np;
      }
    }

    if(Nhalo){
      // send buffer for outgoing halo
      iint *sendBuffer = (iint*) calloc(Nhalo, sizeof(iint));

      // exchange node numbers with neighbors
      meshHaloExchange(mesh,
                     mesh->Np*sizeof(iint),
                     localNums,
                     sendBuffer,
                     localNums+Nlocal);
    }

    preconGatherInfo_t *preconGatherInfoDg =
      (preconGatherInfo_t*) calloc(NpP*mesh->Nelements,
                                 sizeof(preconGatherInfo_t));

    // set local ids
    for(iint n=0;n<mesh->Nelements*NpP;++n)
      preconGatherInfoDg[n].localId = n;

    // numbering of patch interior nodes
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint n=0;n<mesh->Np;++n){
        iint id  = n + e*mesh->Np;
        iint pid = n + e*NpP;

        // all patch interior nodes are local
        preconGatherInfoDg[pid].baseId = localNums[id];
      }
    }

    // add patch boundary nodes
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint f=0;f<mesh->Nfaces;++f){
        // mark halo nodes
        iint rP = mesh->EToP[e*mesh->Nfaces+f];
        iint eP = mesh->EToE[e*mesh->Nfaces+f];
        iint fP = mesh->EToF[e*mesh->Nfaces+f];
        iint bc = mesh->EToB[e*mesh->Nfaces+f];

        for(iint n=0;n<mesh->Nfp;++n){
        iint id = n + f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces;
        iint idP = mesh->vmapP[id];

        // local numbers
        iint pidM = e*NpP + mesh->faceNodes[f*mesh->Nfp+n];
        iint pidP = e*NpP + mesh->Np + f*mesh->Nfp+n;
        preconGatherInfoDg[pidP].baseId = localNums[idP];

        if(rP!=-1){
          preconGatherInfoDg[pidM].haloFlag = 1;
          preconGatherInfoDg[pidP].haloFlag = 1;
        }
        }
      }
    }

    // sort by rank then base index
    qsort(preconGatherInfoDg, NpP*mesh->Nelements, sizeof(preconGatherInfo_t),
        parallelCompareBaseId);

    // do not gather-scatter nodes labelled zero
    int skip = 0;

    while(preconGatherInfoDg[skip].baseId==0 && skip<NpP*mesh->Nelements){
      ++skip;
    }

    // reset local ids
    iint NlocalDg = NpP*mesh->Nelements - skip;
    iint *gatherLocalIdsDg  = (iint*) calloc(NlocalDg, sizeof(iint));
    iint *gatherBaseIdsDg   = (iint*) calloc(NlocalDg, sizeof(iint));
    iint *gatherHaloFlagsDg = (iint*) calloc(NlocalDg, sizeof(iint));
    for(iint n=0;n<NlocalDg;++n){
      gatherLocalIdsDg[n]  = preconGatherInfoDg[n+skip].localId;
      gatherBaseIdsDg[n]   = preconGatherInfoDg[n+skip].baseId;
      gatherHaloFlagsDg[n] = preconGatherInfoDg[n+skip].haloFlag;
    }

    // make preconBaseIds => preconNumbering
    precon->ogsDg = meshParallelGatherScatterSetup(mesh,
                                                 NlocalDg,
                                                 sizeof(dfloat),
                                                 gatherLocalIdsDg,
                                                 gatherBaseIdsDg,
                                                 gatherHaloFlagsDg);

    //correction for full patch
    NpP = mesh->NpP;

    // build degree vector
    iint NtotalDGP = NpP*mesh->Nelements;
    dfloat *invDegree = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
    dfloat *degree    = (dfloat*) calloc(NtotalDGP, sizeof(dfloat));
    precon->o_invDegreeDGP = mesh->device.malloc(NtotalDGP*sizeof(dfloat), invDegree);

    for(iint n=0;n<NtotalDGP;++n)
      degree[n] = 1;

    occa::memory o_deg = mesh->device.malloc(NtotalDGP*sizeof(dfloat), degree);
    meshParallelGatherScatter(mesh, precon->ogsDg, o_deg, o_deg, dfloatString, "add");
    o_deg.copyTo(degree);
    mesh->device.finish();
    o_deg.free();

    for(iint n=0;n<NtotalDGP;++n){ // need to weight inner products{
      if(degree[n] == 0) printf("WARNING!!!!\n");
      invDegree[n] = 1./degree[n];
    }

    precon->o_invDegreeDGP.copyFrom(invDegree);
    free(degree);
    free(invDegree);

    // -------------------------------------------------------------------------------------------
    // load prebuilt transform matrices
    precon->o_oasForwardDg = mesh->device.malloc(NpP*NpP*sizeof(dfloat), mesh->oasForwardDg);
    precon->o_oasBackDg    = mesh->device.malloc(NpP*NpP*sizeof(dfloat), mesh->oasBackDg);

    dfloat *oasForwardDgT = (dfloat*) calloc(NpP*NpP, sizeof(dfloat));
    dfloat *oasBackDgT = (dfloat*) calloc(NpP*NpP, sizeof(dfloat));
    for(iint n=0;n<NpP;++n){
      for(iint m=0;m<NpP;++m){
        oasForwardDgT[n+m*NpP] = mesh->oasForwardDg[m+n*NpP];
        oasBackDgT[n+m*NpP] = mesh->oasBackDg[m+n*NpP];
      }
    }

    precon->o_oasForwardDgT = mesh->device.malloc(NpP*NpP*sizeof(dfloat), oasForwardDgT);
    precon->o_oasBackDgT    = mesh->device.malloc(NpP*NpP*sizeof(dfloat), oasBackDgT);


    /// ---------------------------------------------------------------------------

    // hack estimate for Jacobian scaling

    dfloat *diagInvOp = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
    dfloat *diagInvOpDg = (dfloat*) calloc(NpP*mesh->Nelements, sizeof(dfloat));
    for(iint e=0;e<mesh->Nelements;++e){

      // S = Jabc*(wa*wb*wc*lambda + wb*wc*Da'*wa*Da + wa*wc*Db'*wb*Db + wa*wb*Dc'*wc*Dc)
      // S = Jabc*wa*wb*wc*(lambda*I+1/wa*Da'*wa*Da + 1/wb*Db'*wb*Db + 1/wc*Dc'*wc*Dc)

      dfloat J = mesh->vgeo[e*mesh->Nvgeo + JID];
      dfloat rx = mesh->vgeo[e*mesh->Nvgeo + RXID];
      dfloat sx = mesh->vgeo[e*mesh->Nvgeo + SXID];
      dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
      dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];

      //metric tensor on this element
      dfloat grr = rx*rx+ry*ry;
      dfloat grs = rx*sx+ry*sy;
      dfloat gss = sx*sx+sy*sy;

      //eigenvalues of the metric
      dfloat eigG1 = 0.5*(grr + gss) - 0.5*sqrt((grr-gss)*(grr-gss) + 4*grs*grs);
      dfloat eigG2 = 0.5*(grr + gss) + 0.5*sqrt((grr-gss)*(grr-gss) + 4*grs*grs);

      //hinv scaling?
      dfloat invhMin = 1e9;
      dfloat invhMax = 0;
      for (iint f=0;f<mesh->Nfaces;f++) {
        iint base = e*mesh->Nfaces + f;
        //invh = mymax(mesh->sgeo[base*mesh->Nsgeo+IHID],invh);
        invhMax =  mymax(mesh->sgeo[base*mesh->Nsgeo + SJID],invhMax);
        invhMin =  mymin(mesh->sgeo[base*mesh->Nsgeo + SJID],invhMin);
      }
      dfloat invh = invhMax/invhMin +invhMin/invhMax;

      //TODO Average? Min/Max/Avg Eigenvalue? What works best for the scaling?
      //dfloat Jhinv2 = J*(eigG1+eigG2);
      dfloat Jhinv2 = (2*mesh->N+1)*J*mymax(eigG1,eigG2);
      //dfloat Jhinv2 = J*mymin(eigG1,eigG2);
      //dfloat Jhinv2 = mesh->N*invh;

      for(iint n=0;n<NpP;++n){
        iint pid = n + e*NpP;

        diagInvOpDg[pid] = 1./(J*lambda + Jhinv2*mesh->oasDiagOpDg[n]);
      }
    }

    precon->o_oasDiagInvOpDg =
      mesh->device.malloc(NpP*mesh->Nelements*sizeof(dfloat), diagInvOpDg);

    if (strstr(options,"PATCHSOLVE")||strstr(options,"APPROXPATCH")
      ||strstr(options,"LOCALPATCH")) {
      iint nnz;
      nonZero_t *A;
      dfloat *invAP;
      dfloat *localInvA;
      hgs_t *hgs;

      NpP = mesh->Np*(mesh->Nfaces+1);
      iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

      //initialize the full inverse operators on each 4 element patch
      ellipticBuildPatchesIpdgTri2D(mesh, mesh->Np, NULL, tau, lambda,
                                   BCType, &A, &nnz, &hgs, globalStarts,
                                   &invAP, &localInvA, options);

      if (strstr(options,"PATCHSOLVE")) {
        precon->o_invAP = mesh->device.malloc(mesh->Nelements*NpP*NpP*sizeof(dfloat),invAP);
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
        precon->o_invAP = mesh->device.malloc(mesh->Nelements*NpP*NpP*sizeof(dfloat),invAP);

        //rotated node ids of neighbouring element
        mesh->o_rmapP = mesh->device.malloc(mesh->Np*mesh->Nfaces*sizeof(iint),mesh->rmapP);
        mesh->o_EToF = mesh->device.malloc(mesh->Nelements*mesh->Nfaces*sizeof(iint),mesh->EToF);
      } else if (strstr(options,"LOCALPATCH")) {
        //initialize the full inverse operators on each 4 element patch
        ellipticBuildPatchesIpdgTri2D(mesh, mesh->Np, NULL, tau, lambda,
                                   BCType, &A, &nnz, &hgs, globalStarts,
                                   &invAP, &localInvA, options);

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
    ellipticCoarsePreconditionerSetupTri2D(mesh, precon, tau, lambda, BCType, options);
    occaTimerToc(mesh->device,"CoarsePreconditionerSetup");

  } else if(strstr(options, "FULLALMOND")){
    iint nnz;
    nonZero_t *A;
    hgs_t *hgs;

    iint Nnum = mesh->Np*mesh->Nelements;
    iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));

    if (strstr(options,"IPDG")) {

      MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
      for(iint r=0;r<size;++r)
        globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Np;

      ellipticBuildIpdgTri2D(mesh, tau, lambda, BCType, &A, &nnz,&hgs,globalStarts, options);

      qsort(A, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

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

    precon->parAlmond = parAlmondSetup(mesh,
                                   globalStarts,
                                   nnz,
                                   Rows,
                                   Cols,
                                   Vals,
                                   hgs,
                                   options);

    free(A); free(Rows); free(Cols); free(Vals);
  }
  else if (strstr(options, "BLOCKJACOBI")){

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
  }


  return precon;
}
