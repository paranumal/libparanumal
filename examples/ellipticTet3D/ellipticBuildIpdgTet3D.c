
#include "ellipticTet3D.h"

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgTet3D(mesh3D *mesh, dfloat tau, dfloat lambda, int *BCType, nonZero_t **A,
                            dlong *nnzA, hlong *globalStarts, const char *options){

  int size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  // number of degrees of freedom on this rank
  hlong Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, MPI_COMM_WORLD);
    for(int r=0;r<size;++r)
      globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(size+1, sizeof(hlong));
  dlong Nelements = mesh->Nelements;
  MPI_Allgather(&(mesh->Nelements), 1, MPI_DLONG,
                      rankNelements, 1, MPI_DLONG, MPI_COMM_WORLD);
  //find offsets
  for(int r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) {
      globalIds[e*mesh->Np +n] = n + (e + rankStarts[rankM])*mesh->Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, mesh->Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  dfloat *qmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *qmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;
        for (int i=0;i<mesh->Np;i++)
          MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]  = MSnm;
      }
    }
  }

  // DrT*MS, DsT*MS, DtT*MS
  dfloat *DrTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DsTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DtTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Np;n++) {
      for (int i=0;i<mesh->Nfp;i++) {
        DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DtTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        for (int m=0;m<mesh->Np;m++) {
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
  dlong nnz = 0;

  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // loop over all elements
  #pragma omp parallel for
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    dlong gbase = eM*mesh->Nggeo;
    dfloat Grr = mesh->ggeo[gbase+G00ID];
    dfloat Grs = mesh->ggeo[gbase+G01ID];
    dfloat Grt = mesh->ggeo[gbase+G02ID];
    dfloat Gss = mesh->ggeo[gbase+G11ID];
    dfloat Gst = mesh->ggeo[gbase+G12ID];
    dfloat Gtt = mesh->ggeo[gbase+G22ID];
    dfloat J   = mesh->ggeo[gbase+GWJID];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
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

    dlong vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat drdz = mesh->vgeo[vbase+RZID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat dsdz = mesh->vgeo[vbase+SZID];
    dfloat dtdx = mesh->vgeo[vbase+TXID];
    dfloat dtdy = mesh->vgeo[vbase+TYID];
    dfloat dtdz = mesh->vgeo[vbase+TZID];

    for (int m=0;m<mesh->Np;m++) {
      for (int fM=0;fM<mesh->Nfaces;fM++) {
        // load surface geofactors for this face
        dlong sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat nz = mesh->sgeo[sid+NZID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat hinv = mesh->sgeo[sid+IHID];

        dlong eP = mesh->EToE[eM*mesh->Nfaces+fM];
        if (eP < 0) eP = eM;
        dlong vbaseP = eP*mesh->Nvgeo;
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
        for (int i=0;i<mesh->Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          int idM    = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
          int vidM   = mesh->faceNodes[i+fM*mesh->Nfp];
          int vidP   = (int) (mesh->vmapP[idM]%mesh->Np); // only use this to identify location of positive trace vgeo

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

        dfloat penalty = tau*hinv;
        eP = mesh->EToE[eM*mesh->Nfaces+fM];
        for (int n=0;n<mesh->Np;n++) {
          for (int i=0;i<mesh->Nfp;i++) {
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
            int bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
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

            for (int i=0;i<mesh->Nfp;i++) {
              BM[m+n*mesh->Np] += -0.5*gradqSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
              BM[m+n*mesh->Np] += +0.5*qSgn*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
              BM[m+n*mesh->Np] += -0.5*qSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
            }
          } else {
            for (int i=0;i<mesh->Nfp;i++) {
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmP[i];
              AnmP += +0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                      +0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                      +0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i];
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmP[i];
            }
          }

          if(fabs(AnmP)>tol){
            #pragma omp critical
            {
              // remote info
              (*A)[nnz].row = globalIds[eM*mesh->Np+n];
              (*A)[nnz].col = globalIds[eP*mesh->Np+m];
              (*A)[nnz].val = AnmP;
              (*A)[nnz].ownerRank = rankM;
              ++nnz;
            }
          }
        }
      }
    }

    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
        dfloat Anm = BM[m+n*mesh->Np];

        if(fabs(Anm)>tol){
          #pragma omp critical
          {
            (*A)[nnz].row = globalIds[eM*mesh->Np+n];
            (*A)[nnz].col = globalIds[eM*mesh->Np+m];
            (*A)[nnz].val = Anm;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
  }
  
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  // free up unused storage
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;
  
  if(rankM==0) printf("done.\n");
  
  free(globalIds);

  free(BM);  free(MS);
  free(DrTMS); free(DsTMS); free(DtTMS);

  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}
