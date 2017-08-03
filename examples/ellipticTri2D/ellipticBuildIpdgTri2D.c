
#include "ellipticTri2D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A, iint *nnzA,
                      hgs_t **hgs, iint *globalStarts, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);


  iint Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(iint));
  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));

  if (strstr(options,"PROJECT")) {
    // Create a contiguous numbering system, starting from the element-vertex connectivity
    for (iint n=0;n<Nnum;n++) {
      iint id = mesh->gatherLocalIds[n];
      globalIds[id] = mesh->gatherBaseIds[n];
    }

    // squeeze node numbering
    meshParallelConsecutiveGlobalNumbering(Nnum, globalIds, globalOwners, globalStarts);

    //use the ordering to define a gather+scatter for assembly
    *hgs = meshParallelGatherSetup(mesh, Nnum, globalIds, globalOwners);

  } else {
    // every degree of freedom has its own global id
    MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
      for(iint r=0;r<size;++r)
        globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Np;

    /* so find number of elements on each rank */
    iint *rankNelements = (iint*) calloc(size, sizeof(iint));
    iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
    MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
      rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
    //find offsets
    for(iint r=0;r<size;++r){
      rankStarts[r+1] = rankStarts[r]+rankNelements[r];
    }
    //use the offsets to set a global id
    for (iint e =0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Np;n++) {
        globalIds[e*mesh->Np +n] = n + (e + rankStarts[rankM])*mesh->Np;
        globalOwners[e*mesh->Np +n] = rankM;
      }
    }

    /* do a halo exchange of global node numbers */
    if (mesh->totalHaloPairs) {
      iint *idSendBuffer = (iint *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(iint));
      meshHaloExchange(mesh, mesh->Np*sizeof(iint), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
      free(idSendBuffer);
    }
  }

  iint nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  int forceSymmetry = (strstr(options, "FORCESYMMETRY")) ? 1:0;
  dfloat scale = 1;
  if(forceSymmetry){
    scale = 0.5;
    nnzLocalBound *= 2;
  }

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-12;

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  dfloat *qmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *qmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Nfp;m++) {
  dfloat MSnm = 0;

  for (iint i=0;i<mesh->Np;i++)
    MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];

  MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]  = MSnm;
      }
    }
  }

  // DrT*MS, DsT*MS
  dfloat *DrTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DsTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint i=0;i<mesh->Nfp;i++) {
        DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        for (iint m=0;m<mesh->Np;m++) {
          DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Dr[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Ds[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
        }
      }
    }
  }

  // reset non-zero counter
  int nnz = 0;

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    iint gbase = eM*mesh->Nggeo;
    dfloat Grr = mesh->ggeo[gbase+G00ID];
    dfloat Grs = mesh->ggeo[gbase+G01ID];
    dfloat Gss = mesh->ggeo[gbase+G11ID];
    dfloat J   = mesh->ggeo[gbase+GWJID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
        BM[m+n*mesh->Np]  = J*lambda*mesh->MM[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grr*mesh->Srr[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grs*mesh->Srs[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grs*mesh->Ssr[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gss*mesh->Sss[m+n*mesh->Np];
      }
    }

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];

    for (iint m=0;m<mesh->Np;m++) {
      for (iint fM=0;fM<mesh->Nfaces;fM++) {
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat hinv = mesh->sgeo[sid+IHID];

        iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
        if (eP < 0) eP = eM;
        iint vbaseP = eP*mesh->Nvgeo;
        dfloat drdxP = mesh->vgeo[vbaseP+RXID];
        dfloat drdyP = mesh->vgeo[vbaseP+RYID];
        dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
        dfloat dsdyP = mesh->vgeo[vbaseP+SYID];

        // extract trace nodes
        for (iint i=0;i<mesh->Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          iint idM    = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
          iint vidM   = mesh->faceNodes[i+fM*mesh->Nfp];
          iint vidP   = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo

          qmM[i] =0;
          if (vidM == m) qmM[i] =1;
          qmP[i] =0;
          if (vidP == m) qmP[i] =1;

          ndotgradqmM[i] = (nx*drdx+ny*drdy)*mesh->Dr[m+vidM*mesh->Np]
                          +(nx*dsdx+ny*dsdy)*mesh->Ds[m+vidM*mesh->Np];
          ndotgradqmP[i] = (nx*drdxP+ny*drdyP)*mesh->Dr[m+vidP*mesh->Np]
                          +(nx*dsdxP+ny*dsdyP)*mesh->Ds[m+vidP*mesh->Np];
        }

        dfloat penalty = tau*hinv; // tau*(mesh->N+1)*(mesh->N+1)*hinv;
        eP = mesh->EToE[eM*mesh->Nfaces+fM];


        for (iint n=0;n<mesh->Np;n++) {
          for (iint i=0;i<mesh->Nfp;i++) {
            BM[m+n*mesh->Np] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
            BM[m+n*mesh->Np] += -0.5*sJ*(nx*drdx+ny*drdy)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                -0.5*sJ*(nx*dsdx+ny*dsdy)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
            BM[m+n*mesh->Np] += +0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
          }

          dfloat AnmP = 0;
          if (eP < 0) {
            int qSgn, gradqSgn;
            int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
            iint bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
            if(bcType==1){ // Dirichlet
              qSgn     = -1;
              gradqSgn =  1;
            } else if (bcType==2){ // Neumann
              qSgn     =  1;
              gradqSgn = -1;
            } else { // Neumann for now
              qSgn     =  1;
              gradqSgn = -1;
            }

            for (iint i=0;i<mesh->Nfp;i++) {
              BM[m+n*mesh->Np] += -0.5*gradqSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
              BM[m+n*mesh->Np] +=
                  +0.5*qSgn*sJ*(nx*drdx+ny*drdy)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                  +0.5*qSgn*sJ*(nx*dsdx+ny*dsdy)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
              BM[m+n*mesh->Np] += -0.5*qSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
            }
          } else {
            for (iint i=0;i<mesh->Nfp;i++) {
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmP[i];
              AnmP +=
                  +0.5*sJ*(nx*drdx+ny*drdy)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                  +0.5*sJ*(nx*dsdx+ny*dsdy)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i];
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmP[i];
            }
          }

          if(fabs(AnmP)>tol){
            sendNonZeros[nnz].row = globalIds[eM*mesh->Np + n];
            sendNonZeros[nnz].col = globalIds[eP*mesh->Np + m];
            sendNonZeros[nnz].val = scale*AnmP;
            sendNonZeros[nnz].ownerRank = globalOwners[eM*mesh->Np + n];
            ++nnz;

            if(forceSymmetry){
              sendNonZeros[nnz].row = globalIds[eP*mesh->Np + m];
              sendNonZeros[nnz].col = globalIds[eM*mesh->Np + n];
              sendNonZeros[nnz].val = scale*AnmP;
              sendNonZeros[nnz].ownerRank = globalOwners[eP*mesh->Np + m];
              ++nnz;
            }
          }
        }
      }
    }

    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        dfloat Anm = BM[m+n*mesh->Np];

        if(fabs(Anm)>tol){
          sendNonZeros[nnz].row = globalIds[eM*mesh->Np+n];
          sendNonZeros[nnz].col = globalIds[eM*mesh->Np+m];
          sendNonZeros[nnz].val = scale*Anm;
          sendNonZeros[nnz].ownerRank = globalOwners[eM*mesh->Np + n];
          ++nnz;

          if(forceSymmetry){
            sendNonZeros[nnz].row = globalIds[eM*mesh->Np+m];
            sendNonZeros[nnz].col = globalIds[eM*mesh->Np+n];
            sendNonZeros[nnz].val = scale*Anm;
            sendNonZeros[nnz].ownerRank = globalOwners[eM*mesh->Np + m];
            ++nnz;
          }
        }
      }
    }
  }

  // count how many non-zeros to send to each process
  for(iint n=0;n<nnz;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_IINT, ArecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnzA = 0;
  for(iint r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnzA += ArecvCounts[r]/sizeof(nonZero_t);
  }

  *A = (nonZero_t*) calloc(*nnzA, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
    (*A), ArecvCounts, ArecvOffsets, MPI_CHAR,
    MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnzA, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  nnz = 0;
  for(iint n=1;n<*nnzA;++n){
    if((*A)[n].row == (*A)[nnz].row &&
       (*A)[n].col == (*A)[nnz].col){
      (*A)[nnz].val += (*A)[n].val;
    }
    else{
      ++nnz;
      (*A)[nnz] = (*A)[n];
    }
  }
  *nnzA = nnz+1;

  free(globalIds);
  free(globalOwners);
  free(sendNonZeros);
  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);

  free(BM);  free(MS);
  free(DrTMS); free(DsTMS);

  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}
