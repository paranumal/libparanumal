
#include "ellipticHex3D.h"

typedef struct{

  iint row;
  iint col;
  dfloat val;

}nonZero_t;

void ellipticBuildIpdgHex3D(mesh3D *mesh, dfloat lambda, dfloat *Avals, iint *Arows, iint *Acols, const char *options){

  /* do a halo exchange of local node numbers */

  iint nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;
  nonZero_t *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  for(iint e=0;e<mesh->Nelements;++e){
    
    /* build Dx,Dy,Dz (forget the TP for the moment) */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
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

	for(iint f=0;f<mesh->Nfaces;++f){
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

	    // form negative trace terms in IPDG
	    int idnM = n*mesh->Np+vidM;
	    int idmM = m*mesh->Np+vidM;
	    dfloat dlndxM = drdxM*Br[idn] + dsdxM*Bs[idn] + dtdxM*Bt[idn];
	    dfloat dlndyM = drdyM*Br[idn] + dsdyM*Bs[idn] + dtdyM*Bt[idn];
	    dfloat dlndzM = drdzM*Br[idn] + dsdzM*Bs[idn] + dtdzM*Bt[idn];	  
	    dfloat dlmdxM = drdxM*Br[idm] + dsdxM*Bs[idm] + dtdxM*Bt[idm];
	    dfloat dlmdyM = drdyM*Br[idm] + dsdyM*Bs[idm] + dtdyM*Bt[idm];
	    dfloat dlmdzM = drdzM*Br[idm] + dsdzM*Bs[idm] + dtdzM*Bt[idm];
	    dfloat lnM = B[idn], lmM = B[idm];
	    dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
	    dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM+nz*dlmdzM;
	    Anm += wsJ*lnM*ndotgradlmM; // (ln^-, N.grad lm^-) check sign
	    Anm += wsJ*ndotgradlnM*lmM; // (N.grad ln^-, lm^-) check sign
	    
	  }
	}

	dfloat tol = 1e-8;
	if(fabs(Anm)>tol){
	  // local block
	  A[cnt].row = n + e*mesh->Np + cumStarts[r];
	  A[cnt].col = m + e*mesh->Np + cumStarts[r];
	  A[cnt].val = Anm;
	  ++cnt;
	}
      }
    }
  }
}
