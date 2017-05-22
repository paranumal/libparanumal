
#include "ellipticHex3D.h"

typedef struct{

  iint row;
  iint col;
  dfloat val;

}nonZero_t;

void ellipticBuildIpdgHex3D(mesh3D *mesh, dfloat lambda, nonZero_t **A, iint *nnzA, const char *options){

  /* do a halo exchange of local node numbers */

  iint nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

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
		Br[mode*mesh->Np+node] = mesh->D[ni*mesh->Nq+i]; // check order
	      if(ni==i && nk==k)
		Bs[mode*mesh->Np+node] = mesh->D[nj*mesh->Nq+j]; // check order
	      if(ni==i && nj==j)
		Bt[mode*mesh->Np+node] = mesh->D[nk*mesh->Nq+k]; // check order
	      
	      ++node;
	    }
	  }
	}
	
	++mode;
      }
    }
  }
	
  // reset non-zero counter
  int nnz = 0;
      
  // loop over all elements
  for(iint e=0;e<mesh->Nelements;++e){
    
    /* build Dx,Dy,Dz (forget the TP for the moment) */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
	dfloat Anm = 0;

	// (grad phi_n, grad phi_m)_{D^e}
	for(iint i=0;i<mesh->Np;++i){
	  iint base = e*mesh->Np*mesh->Nvgeo + i;
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
	  Anm += JW*B[idn]*B[idm];
	}

	// loop over all faces in this element
	for(iint f=0;f<mesh->Nfaces;++f){
	  // accumulate flux terms for negative and positive traces
	  dfloat AnmP = 0;
	  for(iint i=0;i<mesh->Nfp;++i){
	    iint vidM = mesh->faceNodes[i+f*mesh->Nfp];

	    // grab vol geofacs at surface nodes
	    iint baseM = e*mesh->Np*mesh->Nvgeo + vidM;
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
	    iint vidP = mesh->vmapP[e*mesh->Nfp*mesh->nNfaces+f*mesh->Nfp+i]%mesh->Np;
	    iint baseP = eP*mesh->Np*mesh->Nvgeo + vidP;
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
	    iint base = mesh->Nsgeo*(e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + i);
	    dfloat nx = mesh->sgeo[base+NXID];
	    dfloat ny = mesh->sgeo[base+NYID];
	    dfloat nz = mesh->sgeo[base+NZID];
	    dfloat wsJ = mesh->sgeo[base+WSJID];
	    dfloat hinv = mesh->sgeo[base+IHID];
	    
	    // form negative trace terms in IPDG
	    int idnM = n*mesh->Np+vidM; // sort this out
	    int idmM = m*mesh->Np+vidM;
	    int idmP = m*mesh->Np+vidP;

	    dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
	    dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
	    dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];	  
	    dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM] + dtdxM*Bt[idmM];
	    dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM] + dtdyM*Bt[idmM];
	    dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM] + dtdzM*Bt[idmM];

	    dfloat dlndxP = drdxP*Br[idnP] + dsdxP*Bs[idnP] + dtdxP*Bt[idnP];
	    dfloat dlndyP = drdyP*Br[idnP] + dsdyP*Bs[idnP] + dtdyP*Bt[idnP];
	    dfloat dlndzP = drdzP*Br[idnP] + dsdzP*Bs[idnP] + dtdzP*Bt[idnP];	  
	    dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP] + dtdxP*Bt[idmP];
	    dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP] + dtdyP*Bt[idmP];
	    dfloat dlmdzP = drdzP*Br[idmP] + dsdzP*Bs[idmP] + dtdzP*Bt[idmP];

	    dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv;
	    
	    dfloat lnM = B[idnM], lmM = B[idmM];
	    dfloat lnP = B[idnP], lmP = B[idmP];
	    
	    dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
	    dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM+nz*dlmdzM;
	    dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP+nz*dlmdzP;

	    Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
	    Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
	    Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

	    // need to check for BC ?
	    AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
	    AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
	    AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
	  }
	  
	  if(fabs(AnmP)>tol){
	    // local block
	    iint rP = mesh->EToP[e*mesh->Nfaces+f]; // check this 
	    (*A)[nnz].row = n + e*mesh->Np + rankStarts[rM];
	    (*A)[nnz].col = m + e*mesh->Np + rankStarts[rP];
	    (*A)[nnz].val = Anm;
	    ++nnz;
	  }
	}
	
	if(fabs(Anm)>tol){
	  // local block
	  (*A)[nnz].row = n + e*mesh->Np + rankStarts[rM];
	  (*A)[nnz].col = m + e*mesh->Np + rankStarts[rM];
	  (*A)[nnz].val = Anm;
	  ++nnz;
	}
      }
    }
  }

  // free up unused storage
  *A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  
  *nnzA = nnz;
  
  free(B);
  free(Br);
  free(Bs);
  free(Bt);
}
