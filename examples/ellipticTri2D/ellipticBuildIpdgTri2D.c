
#include "ellipticTri2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

} nonZero_t;

void ellipticBuildIpdgTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnzA, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  /* find number of elements on each rank */
  iint *rankNelements = (iint*) calloc(size, sizeof(iint));
  iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
    rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
  for(iint r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+(rankNelements[r])*mesh->Np;
  }
  
  /* do a halo exchange of local node numbers */

  iint nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;
  dfloat tau = 2; // hackery

  dfloat *S = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *BP = (dfloat *) calloc(mesh->Np*mesh->Np*mesh->Nfaces,sizeof(dfloat));

  dfloat *q = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dqdr = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dqds = (dfloat *) calloc(mesh->Np,sizeof(dfloat));

  dfloat *qfM = (dfloat *) calloc(mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *qfP = (dfloat *) calloc(mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *qrfM = (dfloat *) calloc(mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *qrfP = (dfloat *) calloc(mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *qsfM = (dfloat *) calloc(mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *qsfP = (dfloat *) calloc(mesh->Np*mesh->Nfp,sizeof(dfloat));

  dfloat *LqM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *LqP = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *DrLqM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *DsLqM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *DrLqP = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *DsLqP = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *LqrM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *LqsM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *LqrP = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *LqsP = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  // reset non-zero counter
  int nnz = 0;
  
  int NfaceNfp = mesh->Nfaces*mesh->Nfp;

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){
    
    iint gbase = eM*mesh->Nggeo;
    dfloat Grr = mesh->ggeo[gbase+G00ID];
    dfloat Grs = mesh->ggeo[gbase+G01ID];
    dfloat Gss = mesh->ggeo[gbase+G11ID];

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J    = mesh->vgeo[vbase+ JID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){ 
        S[m+n*mesh->Np] = J*lambda*mesh->MM[m+n*mesh->Np];
        S[m+n*mesh->Np] += Grr*mesh->Srr[m+n*mesh->Np];
        S[m+n*mesh->Np] += Grs*mesh->Srs[m+n*mesh->Np];
        S[m+n*mesh->Np] += Grs*mesh->Ssr[m+n*mesh->Np];
        S[m+n*mesh->Np] += Gss*mesh->Sss[m+n*mesh->Np];
      }
    }

    // zero out trace contribution
    for (iint n=0;n<mesh->Np;n++) 
      for (iint m=0;m<mesh->Np;m++) 
        BM[m+n*mesh->Np] =0.;

    for (iint fM=0;fM<mesh->Nfaces;fM++) {
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat invJ = mesh->sgeo[sid+IJID];
      dfloat hinv = mesh->sgeo[sid+IHID];


      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) eP = eM;
      iint vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
      dfloat JP    = mesh->vgeo[vbaseP+ JID];

      // zero out trace contribution from neighbour
      for (iint n=0;n<mesh->Np;n++) 
        for (iint m=0;m<mesh->Np;m++) 
          BP[m+n*mesh->Np] =0.;

      /* build traces */
      for(iint m=0;m<mesh->Np;++m) {// m will be the sub-block index for negative and positive trace
        
        //build basis vectors
        for (iint i=0;i<mesh->Np;i++) 
          q[i] = 0;
        
        for (iint i=0;i<mesh->Np;i++) {
          dqdr[i] = mesh->Dr[i*mesh->Np+m];
          dqds[i] = mesh->Ds[i*mesh->Np+m];
        }

        // collect traces
        for (iint i=0;i<mesh->Nfp;i++) {
          iint vidM = mesh->faceNodes[i+fM*mesh->Nfp];

          // double check vol geometric factors are in halo storage of vgeo
          iint idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
          iint vidP    = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo

          qfM[i+m*mesh->Nfp] = q[vidM];
          qfP[i+m*mesh->Nfp] = q[vidP];

          qrfM[i+m*mesh->Nfp] = dqdr[vidM];
          qrfP[i+m*mesh->Nfp] = dqdr[vidP];

          qsfM[i+m*mesh->Nfp] = dqds[vidM];
          qsfP[i+m*mesh->Nfp] = dqds[vidP];
        }
      }

      // LIFT*q, LIFT*dqdr, LIFT*dqds
      for (iint n=0;n<mesh->Np;n++) {
        for (iint m=0;m<mesh->Np;m++) {
          LqM[m+n*mesh->Np] = 0.;
          LqP[m+n*mesh->Np] = 0.;
          LqrM[m+n*mesh->Np] = 0.;
          LqrP[m+n*mesh->Np] = 0.;
          LqsM[m+n*mesh->Np] = 0.;
          LqsP[m+n*mesh->Np] = 0.;  
          for (iint i=0;i<mesh->Nfp;i++) {
            LqM[m+n*mesh->Np]  += mesh->LIFT[i+fM*mesh->Nfp+n*NfaceNfp]*qfM [i+m*mesh->Nfp];
            LqP[m+n*mesh->Np]  += mesh->LIFT[i+fM*mesh->Nfp+n*NfaceNfp]*qfP [i+m*mesh->Nfp];
            LqrM[m+n*mesh->Np] += mesh->LIFT[i+fM*mesh->Nfp+n*NfaceNfp]*qrfM[i+m*mesh->Nfp];
            LqrP[m+n*mesh->Np] += mesh->LIFT[i+fM*mesh->Nfp+n*NfaceNfp]*qrfP[i+m*mesh->Nfp];
            LqsM[m+n*mesh->Np] += mesh->LIFT[i+fM*mesh->Nfp+n*NfaceNfp]*qsfM[i+m*mesh->Nfp];
            LqsP[m+n*mesh->Np] += mesh->LIFT[i+fM*mesh->Nfp+n*NfaceNfp]*qsfP[i+m*mesh->Nfp];
          }
        }
      }

      // DrT*LIFTq, DsT*LIFTq
      for (iint n=0;n<mesh->Np;n++) {
        for (iint m=0;m<mesh->Np;m++) {
          DrLqM[m+n*mesh->Np] = 0.;
          DsLqM[m+n*mesh->Np] = 0.;
          DrLqP[m+n*mesh->Np] = 0.;
          DsLqP[m+n*mesh->Np] = 0.;
          for (iint i=0;i<mesh->Np;i++) {
            DrLqM[m+n*mesh->Np] += mesh->Dr[n+i*mesh->Np]*LqM[m+i*mesh->Np];
            DsLqM[m+n*mesh->Np] += mesh->Ds[n+i*mesh->Np]*LqM[m+i*mesh->Np];
            DrLqP[m+n*mesh->Np] += mesh->Dr[n+i*mesh->Np]*LqP[m+i*mesh->Np];
            DsLqP[m+n*mesh->Np] += mesh->Ds[n+i*mesh->Np]*LqP[m+i*mesh->Np];
          }
        }
      }  

      // collect boundary terms
      dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv;    
      for (iint n=0;n<mesh->Np;n++) {
        for (iint m=0;m<mesh->Np;m++) {
          BM[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*drdx+ny*drdy)*DrLqM[m+n*mesh->Np];
          BM[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*dsdx+ny*dsdy)*DsLqM[m+n*mesh->Np];
          BM[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*drdx+ny*drdy)*LqrM[m+n*mesh->Np];
          BM[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*dsdx+ny*dsdy)*LqsM[m+n*mesh->Np];
          BM[m+n*mesh->Np] += +0.5*sJ*invJ*penalty*LqM[m+n*mesh->Np];
        }
      }
      eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) { //Neumann boundary condition
        for (iint n=0;n<mesh->Np;n++) {
          for (iint m=0;m<mesh->Np;m++) { 
            BM[m+n*mesh->Np] += +0.5*sJ*invJ*(nx*drdx+ny*drdy)*DrLqM[m+n*mesh->Np];
            BM[m+n*mesh->Np] += +0.5*sJ*invJ*(nx*dsdx+ny*dsdy)*DsLqM[m+n*mesh->Np];
            BM[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*drdx+ny*drdy)*LqrM[m+n*mesh->Np];
            BM[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*dsdx+ny*dsdy)*LqsM[m+n*mesh->Np];
            BM[m+n*mesh->Np] += -0.5*sJ*invJ*penalty*LqM[m+n*mesh->Np];
          }
        }
      } else {
        for (iint n=0;n<mesh->Np;n++) {
          for (iint m=0;m<mesh->Np;m++) {
            BP[m+n*mesh->Np] += +0.5*sJ*invJ*(nx*drdxP+ny*drdyP)*DrLqP[m+n*mesh->Np];
            BP[m+n*mesh->Np] += +0.5*sJ*invJ*(nx*dsdxP+ny*dsdyP)*DsLqP[m+n*mesh->Np];
            BP[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*drdxP+ny*drdyP)*LqrP[m+n*mesh->Np];
            BP[m+n*mesh->Np] += -0.5*sJ*invJ*(nx*dsdxP+ny*dsdyP)*LqsP[m+n*mesh->Np];
            BP[m+n*mesh->Np] += -0.5*sJ*invJ*penalty*LqP[m+n*mesh->Np];
          }
        }
      }

      for (iint n=0;n<mesh->Np;n++) {
        for (iint m=0;m<mesh->Np;m++) {
          dfloat AnmP =0.;
          for (iint i=0;i<mesh->Np;i++) 
            AnmP += J*mesh->MM[i+n*mesh->Np]*BP[m+i*mesh->Np];

          if(fabs(AnmP)>tol){
            // remote info
            iint rankP = mesh->EToP[eM*mesh->Nfaces+fM]; 
            if (rankP<0) rankP = rankM;
            (*A)[nnz].row = n + eM*mesh->Np + rankStarts[rankM];
            (*A)[nnz].col = m + eP*mesh->Np + rankStarts[rankP]; // this is wrong
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        } 
      }
    }

    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
        dfloat Anm = S[m+n*mesh->Np];
        for (iint i=0;i<mesh->Np;i++) 
          Anm += J*mesh->MM[i+n*mesh->Np]*BM[m+i*mesh->Np];

        if(fabs(Anm)>tol){
          (*A)[nnz].row = n + eM*mesh->Np + rankStarts[rankM];
          (*A)[nnz].col = m + eM*mesh->Np + rankStarts[rankM]; // this is wrong
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      } 
    }
  }

  // free up unused storage
  *A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  
  *nnzA = nnz;
  

  free(BM);  free(BP);  

  free(q);   free(dqdr);  free(dqds);
  free(qfM);  free(qfP);  free(qrfM);  
  free(qrfP); free(qsfM); free(qsfP);
  free(LqM);  free(LqP);  
  free(DrLqM);   free(DsLqM); free(DrLqP);  free(DsLqP);  
  free(LqrM);  free(LqsM);  free(LqrP);  free(LqsP);  
  free(rankNelements); free(rankStarts);
}
