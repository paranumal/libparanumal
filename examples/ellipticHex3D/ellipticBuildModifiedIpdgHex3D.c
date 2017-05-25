
#include "ellipticHex3D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

} nonZero_t;

extern "C"
{
  void sgeev_(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA, float *WR, float *WI,
	      float *VL, int *LDVL, float *VR, int *LDVR, float *WORK, int *LWORK, int *INFO );
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
	      double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
}

void geev(char *JOBVL, char *JOBVR, int *N, dfloat *A, int *LDA, dfloat *WR, dfloat *WI,
	  dfloat *VL, int *LDVL, dfloat *VR, int *LDVR, dfloat *WORK, int *LWORK, int *INFO ){
  
#if dfloatSize==4
  sgeev_
#else
    dgeev_
#endif
    (JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDA, VR, LDA, WORK, LWORK, INFO);
  
}

void eig(const int Nrows, dfloat *A, dfloat *VL, dfloat *VR, dfloat *WR, dfloat *WI){

  int NB  = 256;
  char JOBVL  = 'V';
  char JOBVR  = 'V';
  int     N = Nrows;
  int   LDA = Nrows;
  int  LWORK  = (NB+2)*N;

  dfloat *WORK = (dfloat*) calloc(LWORK, sizeof(dfloat));
  
  int INFO = -999;

  geev(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
       VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);

  free(WORK);
}


void ellipticBuildModifiedIpdgHex3D(mesh3D *mesh,
				    dfloat lambda, nonZero_t **A, iint *nnzA, dfloat **VR,
				    const char *options){

  // set up change of basis:
  
  // split A into  (DD + (A-DD))*x = b  where DD is the block diagonal of A
  // build VR such that  VR'*DD*VR = LL where LL is the diagonal of eigs of DD
  // [ using DD*VR = VR*LL ]
  // i.e.  (VR*LL*VR' + (A-DD))*x = b
  // becomes  VR*(LL + VR'*(A-DD)*VR)*VR'*x = b
  // so solve  (LL+VR'*(A-DD)*VR)*x~ = VR'*b
  // then compute x = VR*x~
  
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

  iint nnzLocalBound = mesh->Np*mesh->Nelements + mesh->Np*mesh->Np*mesh->Nfaces*mesh->Nelements;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));
  dfloat *VL = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(nonZero_t));
  dfloat *WR = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(nonZero_t));
  dfloat *WI = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(nonZero_t));
  dfloat *VRxC = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(nonZero_t));
  *VR = (dfloat*) calloc(mesh->Np*mesh->Np*mesh->Nelements, sizeof(nonZero_t));
  
  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;
  dfloat tau = 2; // hackery

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  iint mode = 0;
  for(iint nk=0;nk<mesh->N+1;++nk){
    for(iint nj=0;nj<mesh->N+1;++nj){
      for(iint ni=0;ni<mesh->N+1;++ni){

      	iint node = 0;

      	for(iint k=0;k<mesh->N+1;++k){
      	  for(iint j=0;j<mesh->N+1;++j){
      	    for(iint i=0;i<mesh->N+1;++i){

      	      if(nk==k && nj==j && ni==i)
      		      B[mode*mesh->Np+node] = 1;
      	      if(nj==j && nk==k)
      		      Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
      	      if(ni==i && nk==k)
      		      Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
      	      if(ni==i && nj==j)
      		      Bt[mode*mesh->Np+node] = mesh->D[nk+mesh->Nq*k]; 
      	      
      	      ++node;
      	    }
      	  }
      	}
      	
      	++mode;
      }
    }
  }

  dfloat *C  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  
  // loop over all elements
#pragma omp parallel for
  for(iint eM=0;eM<mesh->Nelements;++eM){

    iint cnt = mesh->Nelements*mesh->Np + eM*mesh->Np*mesh->Np*mesh->Nfaces;
    
    /* build Dx,Dy,Dz (forget the TP for the moment) */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
      	dfloat Anm = 0;

      	// (grad phi_n, grad phi_m)_{D^e}
      	for(iint i=0;i<mesh->Np;++i){
      	  iint base = eM*mesh->Np*mesh->Nvgeo + i;
      	  dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
      	  dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
      	  dfloat drdz = mesh->vgeo[base+mesh->Np*RZID];
      	  dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
      	  dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
      	  dfloat dsdz = mesh->vgeo[base+mesh->Np*SZID];
      	  dfloat dtdx = mesh->vgeo[base+mesh->Np*TXID];
      	  dfloat dtdy = mesh->vgeo[base+mesh->Np*TYID];
      	  dfloat dtdz = mesh->vgeo[base+mesh->Np*TZID];
      	  dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
      	  
      	  int idn = n*mesh->Np+i;
      	  int idm = m*mesh->Np+i;
      	  dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
      	  dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
      	  dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];	  
      	  dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm] + dtdx*Bt[idm];
      	  dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm] + dtdy*Bt[idm];
      	  dfloat dlmdz = drdz*Br[idm] + dsdz*Bs[idm] + dtdz*Bt[idm];
      	  Anm += JW*(dlndx*dlmdx+dlndy*dlmdy+dlndz*dlmdz);
      	  Anm += lambda*JW*B[idn]*B[idm];
      	}

      	// loop over all faces in this element
      	for(iint fM=0;fM<mesh->Nfaces;++fM){
      	  // accumulate flux terms for negative traces
      	  for(iint i=0;i<mesh->Nfp;++i){
      	    iint vidM = mesh->faceNodes[i+fM*mesh->Nfp];
	    
      	    // grab vol geofacs at surface nodes
      	    iint baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
      	    dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
      	    dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
      	    dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
      	    dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
      	    dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
      	    dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];
      	    dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
      	    dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
      	    dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];

	    // double check vol geometric factors are in halo storage of vgeo
      	    iint idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
      	    iint vidP    = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo
      	    iint localEP = mesh->vmapP[idM]/mesh->Np;

	    // grab surface geometric factors
      	    iint base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
      	    dfloat nx = mesh->sgeo[base+NXID];
      	    dfloat ny = mesh->sgeo[base+NYID];
      	    dfloat nz = mesh->sgeo[base+NZID];
      	    dfloat wsJ = mesh->sgeo[base+WSJID];
      	    dfloat hinv = mesh->sgeo[base+IHID];
	    
      	    // form negative trace terms in IPDG
      	    int idnM = n*mesh->Np+vidM; 
      	    int idmP = m*mesh->Np+vidP;
      	    int idmM = m*mesh->Np+vidM;

      	    dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
      	    dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
      	    dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
      	    dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
      	    dfloat lnM = B[idnM];

      	    dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM] + dtdxM*Bt[idmM];
      	    dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM] + dtdyM*Bt[idmM];
      	    dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM] + dtdzM*Bt[idmM];
      	    dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM+nz*dlmdzM;
      	    dfloat lmM = B[idmM];
      	    
      	    dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv;

	    Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)
	    
            iint eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            if (eP<0) {
	      Anm += +0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
	      Anm += +0.5*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
	      Anm += -0.5*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-) 
            } 
	  }
	}
	
	C[n*mesh->Np+m] = Anm;
      }
    }
    
    // going to use C = VR'*diag(WR)*VR in C*x=b -> (VR'*diag(WR)*VR)*x = b
    // [ relies on symmetry to pass direct to dgeev
    eig(mesh->Np, C, VL, VR[0]+eM*mesh->Np*mesh->Np, WR, WI);
    
    for(iint n=0;n<mesh->Np;++n){
      // local block
      (*A)[cnt].row = n + eM*mesh->Np + rankStarts[rankM];
      (*A)[cnt].col = n + eM*mesh->Np + rankStarts[rankM];
      (*A)[cnt].val = WR[n];
      (*A)[cnt].ownerRank = rankM;
      ++cnt;
    }
  }
  
#pragma omp parallel for
  for(iint eM=0;eM<mesh->Nelements;++eM){
    
    iint cnt = mesh->Nelements*mesh->Np + eM*mesh->Np*mesh->Np*mesh->Nfaces;
    
    for(iint fM=0;fM<mesh->Nfaces;++fM){
      iint eP    = mesh->EToE[eM*mesh->Nfaces+fM];
      iint rankP = mesh->EToP[eM*mesh->Nfaces+fM]; 
      if (eP>=0){
	
	for(iint n=0;n<mesh->Np;++n){
	  for(iint m=0;m<mesh->Np;++m){ 
	    
	    // accumulate flux terms for positive traces
	    dfloat AnmP = 0;
	    for(iint i=0;i<mesh->Nfp;++i){
	      
	      iint vidM = mesh->faceNodes[i+fM*mesh->Nfp];
	      
	      // grab vol geofacs at surface nodes
	      iint baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
	      dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
	      dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
	      dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
	      dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
	      dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
	      dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];
	      dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
	      dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
	      dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];
	      
	      // double check vol geometric factors are in halo storage of vgeo
	      iint idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
	      iint vidP    = mesh->vmapP[idM]%mesh->Np; 
	      iint localEP = mesh->vmapP[idM]/mesh->Np;
	      
	      iint baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; 
	      dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
	      dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
	      dfloat drdzP = mesh->vgeo[baseP+mesh->Np*RZID];
	      dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
	      dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];
	      dfloat dsdzP = mesh->vgeo[baseP+mesh->Np*SZID];
	      dfloat dtdxP = mesh->vgeo[baseP+mesh->Np*TXID];
	      dfloat dtdyP = mesh->vgeo[baseP+mesh->Np*TYID];
	      dfloat dtdzP = mesh->vgeo[baseP+mesh->Np*TZID];

	      // grab surface geometric factors
	      iint base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
	      dfloat nx = mesh->sgeo[base+NXID];
	      dfloat ny = mesh->sgeo[base+NYID];
	      dfloat nz = mesh->sgeo[base+NZID];
	      dfloat wsJ = mesh->sgeo[base+WSJID];
	      dfloat hinv = mesh->sgeo[base+IHID];

	      // form negative trace terms in IPDG
	      int idnM = n*mesh->Np+vidM; 
	      int idmM = m*mesh->Np+vidM;
	      int idmP = m*mesh->Np+vidP;
	      
	      dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
	      dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
	      dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
	      dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
	      dfloat lnM = B[idnM];
	    
	      dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP] + dtdxP*Bt[idmP];
	      dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP] + dtdyP*Bt[idmP];
	      dfloat dlmdzP = drdzP*Br[idmP] + dsdzP*Bs[idmP] + dtdzP*Bt[idmP];
	      dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP+nz*dlmdzP;
	      dfloat lmP = B[idmP];
	      
	      dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv;	   
	      
	      AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
	      AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
	      AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
	    
	    C[n*mesh->Np+m] = AnmP;
	  }
	}
	
	// now pre and post multiply by VR_{eM} and VR_{eP}'
	for(iint n=0;n<mesh->Np;++n){
	  for(iint m=0;m<mesh->Np;++m){
	    dfloat res = 0;
	    for(iint i=0;i<mesh->Np;++i){ // column major VR
	      res += VR[0][eM*mesh->Np*mesh->Np+i+n*mesh->Np]*C[i*mesh->Np+m];
	    }
	    VRxC[n*mesh->Np+m] = res;
	  }
	}
	
	for(iint n=0;n<mesh->Np;++n){
	  for(iint m=0;m<mesh->Np;++m){
	    dfloat res = 0;
	    for(iint i=0;i<mesh->Np;++i){ // column major VR
	      res += VRxC[n*mesh->Np+i]*VR[0][eP*mesh->Np*mesh->Np+n+i*mesh->Np];
	    }
	    
	    (*A)[cnt].row = n + eM*mesh->Np + rankStarts[rankM];
	    (*A)[cnt].col = m + eP*mesh->Np + rankStarts[rankP]; // this is wrong
	    (*A)[cnt].val = res;
	    (*A)[cnt].ownerRank = rankM;
	    ++cnt;
	  }
	}
      }
    }
  }

  // non-zero counter
  int nnz = 0;
  for(iint n=0;n<mesh->Np*mesh->Np*(mesh->Nfaces+1)*mesh->Nelements;++n)
    if(fabs((*A)[n].val)>tol)
      (*A)[nnz++] = (*A)[n];
  
  // free up unused storage
  *A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  
  *nnzA = nnz;
  
  free(B);  free(Br); free(Bs); free(Bt);
  free(C); free(VL); free(WR); free(WI); free(VRxC);
  free(rankNelements); free(rankStarts);
}
