#include "ellipticTri2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

// compare on global indices
int parallelCompareRowColumn(const void *a, const void *b){

  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;

}

void ellipticBuildCoarseIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A, iint *nnzA,
                              hgs_t **hgs, iint *globalStarts, const char *options);

void ellipticBuildCoarseContinuousTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnz,
                              hgs_t **hgs, iint *globalStarts, const char* options);

void ellipticCoarsePreconditionerSetupTri2D(mesh_t *mesh, precon_t *precon, dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **V1, nonZero_t **A, iint *nnzA,
                                   hgs_t **hgs, iint *globalStarts, const char *options){

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint Nnum = mesh->Nverts*(mesh->Nelements+mesh->totalHaloPairs);

  // 2. Build coarse grid element basis functions
  *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vr1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vs1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));

  for(iint n=0;n<mesh->Np;++n){
    dfloat rn = mesh->r[n];
    dfloat sn = mesh->s[n];

    (*V1)[0*mesh->Np+n] = -0.5*(rn+sn);
    (*V1)[1*mesh->Np+n] = +0.5*(1.+rn);
    (*V1)[2*mesh->Np+n] = +0.5*(1.+sn);

    Vr1[0*mesh->Np+n] = 0.5*(-1);
    Vr1[1*mesh->Np+n] = 0.5*(+1);
    Vr1[2*mesh->Np+n] = 0;

    Vs1[0*mesh->Np+n] = 0.5*(-1);
    Vs1[1*mesh->Np+n] = 0;
    Vs1[2*mesh->Np+n] = 0.5*(+1);
  }

  //precon->o_Vr1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vr1);
  //precon->o_Vs1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vs1);

  //build coarse grid A
  if (strstr(options,"IPDG")) {
    MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
    for(iint r=0;r<size;++r)
      globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Nverts;

    ellipticBuildCoarseIpdgTri2D(mesh,tau,lambda,BCType,A,nnzA,hgs,globalStarts,options);

    qsort(*A, *nnzA, sizeof(nonZero_t), parallelCompareRowColumn);

  } else if (strstr(options,"CONTINUOUS")) {

    ellipticBuildCoarseContinuousTri2D(mesh,lambda,A,nnzA,hgs,globalStarts,options);
  }
}


void ellipticBuildCoarseContinuousTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnz,
                              hgs_t **hgs, iint *globalStarts, const char* options) {
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ------------------------------------------------------------------------------------
  // 1. Create a contiguous numbering system, starting from the element-vertex connectivity
  iint Nnum = mesh->Nverts*mesh->Nelements;
  iint *globalNumbering = (iint*) calloc(Nnum, sizeof(iint));
  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));

  // use original vertex numbering
  memcpy(globalNumbering, mesh->EToV, Nnum*sizeof(iint));

  // squeeze numbering
  meshParallelConsecutiveGlobalNumbering(Nnum, globalNumbering, globalOwners, globalStarts);

  //use the ordering to define a gather+scatter for assembly
  *hgs = meshParallelGatherSetup(mesh, Nnum, globalNumbering, globalOwners);

  // ------------------------------------------------------------------------------------
  // 2. Build non-zeros of stiffness matrix (unassembled)
  iint nnzLocal = mesh->Nverts*mesh->Nverts*mesh->Nelements;
  iint   *rowsA = (iint*) calloc(nnzLocal, sizeof(iint));
  iint   *colsA = (iint*) calloc(nnzLocal, sizeof(iint));
  dfloat *valsA = (dfloat*) calloc(nnzLocal, sizeof(dfloat));

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));

  dfloat *cV1  = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVr1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVs1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  for(iint n=0;n<mesh->cubNp;++n){
    dfloat rn = mesh->cubr[n];
    dfloat sn = mesh->cubs[n];

    cV1[0*mesh->cubNp+n] = -0.5*(rn+sn);
    cV1[1*mesh->cubNp+n] = +0.5*(1.+rn);
    cV1[2*mesh->cubNp+n] = +0.5*(1.+sn);

    cVr1[0*mesh->cubNp+n] = 0.5*(-1);
    cVr1[1*mesh->cubNp+n] = 0.5*(+1);
    cVr1[2*mesh->cubNp+n] = 0;

    cVs1[0*mesh->cubNp+n] = 0.5*(-1);
    cVs1[1*mesh->cubNp+n] = 0;
    cVs1[2*mesh->cubNp+n] = 0.5*(+1);
  }

  dfloat Srr1[3][3], Srs1[3][3], Ssr1[3][3], Sss1[3][3], MM1[3][3];
  for(iint n=0;n<mesh->Nverts;++n){
    for(iint m=0;m<mesh->Nverts;++m){
      Srr1[n][m] = 0;
      Srs1[n][m] = 0;
      Ssr1[n][m] = 0;
      Sss1[n][m] = 0;
      MM1[n][m] = 0;

      for(iint i=0;i<mesh->cubNp;++i){
        iint idn = n*mesh->cubNp+i;
        iint idm = m*mesh->cubNp+i;
        dfloat cw = mesh->cubw[i];
        Srr1[n][m] += cw*(cVr1[idn]*cVr1[idm]);
        Srs1[n][m] += cw*(cVr1[idn]*cVs1[idm]);
        Ssr1[n][m] += cw*(cVs1[idn]*cVr1[idm]);
        Sss1[n][m] += cw*(cVs1[idn]*cVs1[idm]);
        MM1[n][m] += cw*(cV1[idn]*cV1[idm]);
      }
    }
  }

  printf("Building coarse matrix system\n");
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      for(iint m=0;m<mesh->Nverts;++m){
        dfloat Snm = 0;

        dfloat rx = mesh->vgeo[e*mesh->Nvgeo + RXID];
        dfloat sx = mesh->vgeo[e*mesh->Nvgeo + SXID];
        dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
        dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];
        dfloat J  = mesh->vgeo[e*mesh->Nvgeo +  JID];

        Snm  = J*(rx*rx+ry*ry)*Srr1[n][m];
        Snm += J*(rx*sx+ry*sy)*Srs1[n][m];
        Snm += J*(sx*rx+sy*ry)*Ssr1[n][m];
        Snm += J*(sx*sx+sy*sy)*Sss1[n][m];
        Snm += J*lambda*MM1[n][m];

        dfloat nonZeroThreshold = 1e-7;
        if(fabs(Snm)>nonZeroThreshold) {
          valsA[cnt] = Snm;
          rowsA[cnt] = e*mesh->Nverts+n;
          colsA[cnt] = e*mesh->Nverts+m;

          // pack non-zero
          sendNonZeros[cnt].val = Snm;
          sendNonZeros[cnt].row = globalNumbering[e*mesh->Nverts+n];
          sendNonZeros[cnt].col = globalNumbering[e*mesh->Nverts+m];
          sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Nverts+n];
          ++cnt;
        }
      }
    }
  }

  // count how many non-zeros to send to each process
  for(iint n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_IINT, ArecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnz = 0;
  for(iint r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnz += ArecvCounts[r]/sizeof(nonZero_t);
  }

  *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
    (*A), ArecvCounts, ArecvOffsets, MPI_CHAR,
    MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(iint n=1;n<*nnz;++n){
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col){
      (*A)[cnt].val += (*A)[n].val;
    }
    else{
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  *nnz = cnt+1;
  printf("Done building coarse matrix system\n");

  free(globalNumbering); free(globalOwners);
  free(cVr1); free(cVs1); free(cV1);
  free(rowsA); free(colsA); free(valsA);
  free(sendNonZeros);

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);
}

void ellipticBuildCoarseIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A, iint *nnzA,
                      hgs_t **hgs, iint *globalStarts, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);


  iint Nnum = mesh->Nverts*mesh->Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts,sizeof(iint));
  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));

  if (strstr(options,"PROJECT")||strstr(options,"PRECONC0")) {
    // Create a contiguous numbering system, starting from the element-vertex connectivity
    memcpy(globalIds, mesh->EToV, Nnum*sizeof(iint));

    // squeeze node numbering
    meshParallelConsecutiveGlobalNumbering(Nnum, globalIds, globalOwners, globalStarts);

    //use the ordering to define a gather+scatter for assembly
    *hgs = meshParallelGatherSetup(mesh, Nnum, globalIds, globalOwners);

  } else {
    // every degree of freedom has its own global id
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
      for (int n=0;n<mesh->Nverts;n++) {
        globalIds[e*mesh->Nverts +n] = n + (e + rankStarts[rankM])*mesh->Nverts;
        globalOwners[e*mesh->Nverts +n] = rankM;
      }
    }

    /* do a halo exchange of global node numbers */
    if (mesh->totalHaloPairs) {
      iint *idSendBuffer = (iint *) calloc(mesh->Nverts*mesh->totalHaloPairs,sizeof(iint));
      meshHaloExchange(mesh, mesh->Nverts*sizeof(iint), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Nverts);
      free(idSendBuffer);
    }
  }

  // Build coarse basis restriction
  dfloat *V  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    dfloat rn = mesh->r[n];
    dfloat sn = mesh->s[n];

    V[0*mesh->Np+n] = -0.5*(rn+sn);
    V[1*mesh->Np+n] = +0.5*(1.+rn);
    V[2*mesh->Np+n] = +0.5*(1.+sn);
  }

  iint nnzLocalBound = mesh->Nverts*mesh->Nverts*(1+mesh->Nfaces)*mesh->Nelements;

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-12;

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *BP = (dfloat *) calloc(mesh->Np*mesh->Np*mesh->Nfaces,sizeof(dfloat));

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

    //zero out BP
    for (iint fM=0;fM<mesh->Nfaces;fM++)
      for (iint n=0;n<mesh->Np;n++)
        for (iint m=0;m<mesh->Np;m++)
          BP[m+n*mesh->Np+fM*mesh->Np*mesh->Np] = 0;

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
          BP[m+n*mesh->Np+fM*mesh->Np*mesh->Np] = AnmP;
        }
      }
    }

    //transform to coarse operator Ac = VAV^T
    for (int fM=0;fM<mesh->Nfaces;fM++) {
      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) eP = eM;

      for (int j=0;j<mesh->Nverts;j++) {
        for (int i=0;i<mesh->Nverts;i++) {
          dfloat AnmP = 0;
          for (iint m=0;m<mesh->Np;m++) {
            for (iint n=0;n<mesh->Np;n++) {
              AnmP += V[n+i*mesh->Np]*BP[m+n*mesh->Np+fM*mesh->Np*mesh->Np]*V[m+j*mesh->Np];
            }
          }
          if(fabs(AnmP)>tol){
            sendNonZeros[nnz].row = globalIds[eM*mesh->Nverts + i];
            sendNonZeros[nnz].col = globalIds[eP*mesh->Nverts + j];
            sendNonZeros[nnz].val = AnmP;
            sendNonZeros[nnz].ownerRank = globalOwners[eM*mesh->Nverts + i];
            ++nnz;
          }
        }
      }
    }

    for (int j=0;j<mesh->Nverts;j++) {
      for (int i=0;i<mesh->Nverts;i++) {
        dfloat Anm = 0;
        for (iint m=0;m<mesh->Np;m++) {
          for (iint n=0;n<mesh->Np;n++) {
            Anm += V[n+i*mesh->Np]*BM[m+n*mesh->Np]*V[m+j*mesh->Np];
          }
        }

        if(fabs(Anm)>tol){
          sendNonZeros[nnz].row = globalIds[eM*mesh->Nverts+i];
          sendNonZeros[nnz].col = globalIds[eM*mesh->Nverts+j];
          sendNonZeros[nnz].val = Anm;
          sendNonZeros[nnz].ownerRank = globalOwners[eM*mesh->Nverts + i];
          ++nnz;
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

  free(V);
  free(BP); free(BM);  free(MS);
  free(DrTMS); free(DsTMS);

  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}