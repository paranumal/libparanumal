
#include "ellipticTet3D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgTet3D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A,
                            iint *nnzA, iint *globalStarts, const char *options){

  iint Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(iint));

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
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    iint *idSendBuffer = (iint *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(iint));
    meshHaloExchange(mesh, mesh->Np*sizeof(iint), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
    free(idSendBuffer);
  }

  iint nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

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

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));

  // reset non-zero counter
  int nnz = 0;

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

        dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv;
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
            for (iint i=0;i<mesh->Nfp;i++) {
              BM[m+n*mesh->Np] += +0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
              BM[m+n*mesh->Np] += +0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
              BM[m+n*mesh->Np] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
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

          if(fabs(AnmP)>tol){
            // remote info
            iint rankP = mesh->EToP[eM*mesh->Nfaces+fM];
            if (rankP<0) rankP = rankM;
            (*A)[nnz].row = globalIds[eM*mesh->Np+n];
            (*A)[nnz].col = globalIds[eP*mesh->Np+m];
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }

    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        dfloat Anm = BM[m+n*mesh->Np];

        if(fabs(Anm)>tol){
          (*A)[nnz].row = globalIds[eM*mesh->Np+n];
          (*A)[nnz].col = globalIds[eM*mesh->Np+m];
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }
  printf("nnz = %lu\n", nnz);
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  // free up unused storage
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  free(globalIds);

  free(BM);  free(MS);
  free(DrTMS); free(DsTMS); free(DtTMS);

  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}
