#include "ellipticTet3D.h"

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

void ellipticBuildCoarseIpdgTet3D(mesh3D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A,
                              iint *nnzA, iint *globalStarts, const char *options);

void ellipticBuildCoarseContinuousTet3D(mesh3D *mesh, dfloat lambda, nonZero_t **A, iint *nnz,
                              hgs_t **hgs, iint *globalStarts, const char* options);

void ellipticCoarsePreconditionerSetupTet3D(mesh_t *mesh, precon_t *precon, dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **V1, nonZero_t **A, iint *nnzA,
                                   hgs_t **hgs, iint *globalStarts, const char *options){

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint Nnum = mesh->Nverts*(mesh->Nelements+mesh->totalHaloPairs);

  *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    dfloat rn = mesh->r[n];
    dfloat sn = mesh->s[n];
    dfloat tn = mesh->t[n];

    (*V1)[0*mesh->Np+n] = -0.5*(1.0+rn+sn+tn);
    (*V1)[1*mesh->Np+n] = +0.5*(1.+rn);
    (*V1)[2*mesh->Np+n] = +0.5*(1.+sn);
    (*V1)[3*mesh->Np+n] = +0.5*(1.+tn);
  }

  //build coarse grid A
  if (strstr(options,"IPDG")) {
    MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT, globalStarts+1, 1, MPI_IINT, MPI_COMM_WORLD);
    for(iint r=0;r<size;++r)
      globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Nverts;

    ellipticBuildCoarseIpdgTet3D(mesh,tau,lambda,BCType,A,nnzA,globalStarts,options);

  } else if (strstr(options,"CONTINUOUS")) {
    ellipticBuildCoarseContinuousTet3D(mesh,lambda,A,nnzA,hgs,globalStarts,options);
  }
}

void ellipticBuildCoarseContinuousTet3D(mesh3D *mesh, dfloat lambda, nonZero_t **A, iint *nnz,
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
  hgs_t *hgs = meshParallelGatherSetup(mesh, Nnum, globalNumbering, globalOwners);

  // ------------------------------------------------------------------------------------
  // 2. Build non-zeros of stiffness matrix (unassembled)
  iint nnz = mesh->Nverts*mesh->Nverts*mesh->Nelements;
  iint   *rowsA = (iint*) calloc(nnz, sizeof(iint));
  iint   *colsA = (iint*) calloc(nnz, sizeof(iint));
  dfloat *valsA = (dfloat*) calloc(nnz, sizeof(dfloat));

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnz, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));

  iint cnt = 0;

  dfloat *cV1  = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVr1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVs1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));
  dfloat *cVt1 = (dfloat*) calloc(mesh->cubNp*mesh->Nverts, sizeof(dfloat));

  for(iint n=0;n<mesh->cubNp;++n){
    dfloat rn = mesh->cubr[n];
    dfloat sn = mesh->cubs[n];
    dfloat tn = mesh->cubt[n];

    cV1[0*mesh->cubNp+n] = -0.5*(1.0+rn+sn+tn);
    cV1[1*mesh->cubNp+n] = +0.5*(1.+rn);
    cV1[2*mesh->cubNp+n] = +0.5*(1.+sn);
    cV1[3*mesh->cubNp+n] = +0.5*(1.+tn);

    cVr1[0*mesh->cubNp+n] = 0.5*(-1);
    cVr1[1*mesh->cubNp+n] = 0.5*(+1);
    cVr1[2*mesh->cubNp+n] = 0;
    cVr1[3*mesh->cubNp+n] = 0;

    cVs1[0*mesh->cubNp+n] = 0.5*(-1);
    cVs1[1*mesh->cubNp+n] = 0;
    cVs1[2*mesh->cubNp+n] = 0.5*(+1);
    cVs1[3*mesh->cubNp+n] = 0;

    cVt1[0*mesh->cubNp+n] = 0.5*(-1);
    cVt1[1*mesh->cubNp+n] = 0;
    cVt1[2*mesh->cubNp+n] = 0;
    cVt1[3*mesh->cubNp+n] = 0.5*(+1);
  }

  dfloat Srr1[4][4], Srs1[4][4], Srt1[4][4];
  dfloat Ssr1[4][4], Sss1[4][4], Sst1[4][4];
  dfloat Str1[4][4], Sts1[4][4], Stt1[4][4];
  dfloat MM1[4][4];
  for(iint n=0;n<mesh->Nverts;++n){
    for(iint m=0;m<mesh->Nverts;++m){
      Srr1[n][m] = 0;
      Srs1[n][m] = 0;
      Srt1[n][m] = 0;
      Ssr1[n][m] = 0;
      Sss1[n][m] = 0;
      Sst1[n][m] = 0;
      Str1[n][m] = 0;
      Sts1[n][m] = 0;
      Stt1[n][m] = 0;
      MM1[n][m] = 0;

      for(iint i=0;i<mesh->cubNp;++i){
      	iint idn = n*mesh->cubNp+i;
      	iint idm = m*mesh->cubNp+i;
      	dfloat cw = mesh->cubw[i];
      	Srr1[n][m] += cw*(cVr1[idn]*cVr1[idm]);
      	Srs1[n][m] += cw*(cVr1[idn]*cVs1[idm]);
        Srt1[n][m] += cw*(cVr1[idn]*cVt1[idm]);
      	Ssr1[n][m] += cw*(cVs1[idn]*cVr1[idm]);
        Sss1[n][m] += cw*(cVs1[idn]*cVs1[idm]);
        Sst1[n][m] += cw*(cVs1[idn]*cVt1[idm]);
        Str1[n][m] += cw*(cVt1[idn]*cVr1[idm]);
        Sts1[n][m] += cw*(cVt1[idn]*cVs1[idm]);
        Stt1[n][m] += cw*(cVt1[idn]*cVt1[idm]);
      	MM1[n][m] += cw*(cV1[idn]*cV1[idm]);
      }
    }
  }


  printf("Building coarse matrix system\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      for(iint m=0;m<mesh->Nverts;++m){
      	dfloat Snm = 0;

      	dfloat rx = mesh->vgeo[e*mesh->Nvgeo + RXID];
      	dfloat sx = mesh->vgeo[e*mesh->Nvgeo + SXID];
        dfloat tx = mesh->vgeo[e*mesh->Nvgeo + TXID];
      	dfloat ry = mesh->vgeo[e*mesh->Nvgeo + RYID];
      	dfloat sy = mesh->vgeo[e*mesh->Nvgeo + SYID];
        dfloat ty = mesh->vgeo[e*mesh->Nvgeo + TYID];
        dfloat rz = mesh->vgeo[e*mesh->Nvgeo + RZID];
        dfloat sz = mesh->vgeo[e*mesh->Nvgeo + SZID];
        dfloat tz = mesh->vgeo[e*mesh->Nvgeo + TZID];
      	dfloat J  = mesh->vgeo[e*mesh->Nvgeo +  JID];

      	Snm  = J*(rx*rx+ry*ry+rz*rz)*Srr1[n][m];
      	Snm += J*(rx*sx+ry*sy+rz*sz)*Srs1[n][m];
        Snm += J*(rx*tx+ry*ty+rz*tz)*Srt1[n][m];
        Snm += J*(sx*rx+sy*ry+sz*rz)*Ssr1[n][m];
      	Snm += J*(sx*sx+sy*sy+sz*sz)*Sss1[n][m];
        Snm += J*(sx*tx+sy*ty+sz*tz)*Sst1[n][m];
        Snm += J*(tx*rx+ty*ry+tz*rz)*Str1[n][m];
        Snm += J*(tx*sx+ty*sy+tz*sz)*Sts1[n][m];
        Snm += J*(tx*tx+ty*ty+tz*tz)*Stt1[n][m];
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

void ellipticBuildCoarseIpdgTet3D(mesh3D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A,
                                  iint *nnzA, iint *globalStarts, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);


  iint Nnum = mesh->Nverts*mesh->Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts,sizeof(iint));

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
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    iint *idSendBuffer = (iint *) calloc(mesh->Nverts*mesh->totalHaloPairs,sizeof(iint));
    meshHaloExchange(mesh, mesh->Nverts*sizeof(iint), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Nverts);
    free(idSendBuffer);
  }

  // Build coarse basis restriction
  dfloat *V  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  for(iint n=0;n<mesh->Np;++n){
    dfloat rn = mesh->r[n];
    dfloat sn = mesh->s[n];
    dfloat tn = mesh->t[n];

    V[0*mesh->Np+n] = -0.5*(rn+sn+tn+1.);
    V[1*mesh->Np+n] = +0.5*(1.+rn);
    V[2*mesh->Np+n] = +0.5*(1.+sn);
    V[3*mesh->Np+n] = +0.5*(1.+tn);
  }

  iint nnzLocalBound = mesh->Nverts*mesh->Nverts*(1+mesh->Nfaces)*mesh->Nelements;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

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

  // DrT*MS, DsT*MS, DtT*MS
  dfloat *DrTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DsTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DtTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint i=0;i<mesh->Nfp;i++) {
        DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DtTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        for (iint m=0;m<mesh->Np;m++) {
          DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Dr[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Ds[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DtTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Dt[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
        }
      }
    }
  }

  // reset non-zero counter
  iint nnz = 0;

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    iint gbase = eM*mesh->Nggeo;
    dfloat Grr = mesh->ggeo[gbase+G00ID];
    dfloat Grs = mesh->ggeo[gbase+G01ID];
    dfloat Grt = mesh->ggeo[gbase+G02ID];
    dfloat Gss = mesh->ggeo[gbase+G11ID];
    dfloat Gst = mesh->ggeo[gbase+G12ID];
    dfloat Gtt = mesh->ggeo[gbase+G22ID];
    dfloat J   = mesh->ggeo[gbase+GWJID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
        BM[m+n*mesh->Np]  = J*lambda*mesh->MM[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grr*mesh->Srr[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grs*mesh->Srs[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grt*mesh->Srt[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grs*mesh->Ssr[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gss*mesh->Sss[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gst*mesh->Sst[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grt*mesh->Str[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gst*mesh->Sts[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gtt*mesh->Stt[m+n*mesh->Np];
      }
    }

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat drdz = mesh->vgeo[vbase+RZID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat dsdz = mesh->vgeo[vbase+SZID];
    dfloat dtdx = mesh->vgeo[vbase+TXID];
    dfloat dtdy = mesh->vgeo[vbase+TYID];
    dfloat dtdz = mesh->vgeo[vbase+TZID];

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
        dfloat nz = mesh->sgeo[sid+NZID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat hinv = mesh->sgeo[sid+IHID];

        iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
        if (eP < 0) eP = eM;
        iint vbaseP = eP*mesh->Nvgeo;
        dfloat drdxP = mesh->vgeo[vbaseP+RXID];
        dfloat drdyP = mesh->vgeo[vbaseP+RYID];
        dfloat drdzP = mesh->vgeo[vbaseP+RZID];
        dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
        dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
        dfloat dsdzP = mesh->vgeo[vbaseP+SZID];
        dfloat dtdxP = mesh->vgeo[vbaseP+TXID];
        dfloat dtdyP = mesh->vgeo[vbaseP+TYID];
        dfloat dtdzP = mesh->vgeo[vbaseP+TZID];

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

          ndotgradqmM[i] = (nx*drdx+ny*drdy+nz*drdz)*mesh->Dr[m+vidM*mesh->Np]
                          +(nx*dsdx+ny*dsdy+nz*dsdz)*mesh->Ds[m+vidM*mesh->Np]
                          +(nx*dtdx+ny*dtdy+nz*dtdz)*mesh->Dt[m+vidM*mesh->Np];
          ndotgradqmP[i] = (nx*drdxP+ny*drdyP+nz*drdzP)*mesh->Dr[m+vidP*mesh->Np]
                          +(nx*dsdxP+ny*dsdyP+nz*dsdzP)*mesh->Ds[m+vidP*mesh->Np]
                          +(nx*dtdxP+ny*dtdyP+nz*dtdzP)*mesh->Dt[m+vidP*mesh->Np];
        }

        dfloat penalty = tau*hinv; // tau*(mesh->N+1)*(mesh->N+1)*hinv;
        eP = mesh->EToE[eM*mesh->Nfaces+fM];
        for (iint n=0;n<mesh->Np;n++) {
          for (iint i=0;i<mesh->Nfp;i++) {
            BM[m+n*mesh->Np] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
            BM[m+n*mesh->Np] += -0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                -0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                -0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
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
              BM[m+n*mesh->Np] += +0.5*gradqSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
              BM[m+n*mesh->Np] += +0.5*qSgn*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
              BM[m+n*mesh->Np] += -0.5*qSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
            }
          } else {
            for (iint i=0;i<mesh->Nfp;i++) {
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmP[i];
              AnmP += +0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                      +0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                      +0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i];
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
            (*A)[nnz].row = globalIds[eM*mesh->Nverts + i];
            (*A)[nnz].col = globalIds[eP*mesh->Nverts + j];
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = rankM;
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
          (*A)[nnz].row = globalIds[eM*mesh->Nverts+i];
          (*A)[nnz].col = globalIds[eM*mesh->Nverts+j];
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnzA, sizeof(nonZero_t), parallelCompareRowColumn);

  free(globalIds);

  free(V);
  free(BP); free(BM);  free(MS);
  free(DrTMS); free(DsTMS); free(DtTMS);

  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}