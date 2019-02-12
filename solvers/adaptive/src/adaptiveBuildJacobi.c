/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "adaptive.h"

void BuildLocalIpdgDiagHex3D (adaptive_t* adaptive, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *A);

void BuildLocalContinuousDiagHex3D (adaptive_t* adaptive, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dfloat *A);

void adaptiveBuildJacobi(adaptive_t* adaptive, dfloat lambda, dfloat **invDiagA){

  mesh_t *mesh = adaptive->mesh;
  setupAide options = adaptive->options;

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nfp;n++) {
      int fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (int m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }

  // build some monolithic basis arrays (for quads and hexes)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));


  if (adaptive->elementType==HEXAHEDRA) {
    int mode = 0;
    for(int nk=0;nk<mesh->N+1;++nk){
      for(int nj=0;nj<mesh->N+1;++nj){
        for(int ni=0;ni<mesh->N+1;++ni){

          int node = 0;

          for(int k=0;k<mesh->N+1;++k){
            for(int j=0;j<mesh->N+1;++j){
              for(int i=0;i<mesh->N+1;++i){

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
  }



  dlong diagNnum = mesh->Np*mesh->Nelements;

  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  if(mesh->rank==0) printf("Building diagonal...");fflush(stdout);

  if (options.compareArgs("DISCRETIZATION","IPDG")) {
#pragma omp parallel for 
    for(dlong eM=0;eM<mesh->Nelements;++eM)
      BuildLocalIpdgDiagHex3D(adaptive, mesh, lambda, MS, B, Br, Bs, Bt, eM, diagA + eM*mesh->Np); 
  } else if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
#pragma omp parallel for 
    for(dlong eM=0;eM<mesh->Nelements;++eM)
      BuildLocalContinuousDiagHex3D(adaptive, mesh, lambda, eM, B, Br, Bs, Bt, diagA + eM*mesh->Np);
  }

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) 
    ogsGatherScatter(diagA, ogsDfloat, ogsAdd, adaptive->ogs);
    
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(mesh->rank==0) printf("done.\n");

  free(diagA);
  free(MS);
  free(B); free(Br); free(Bs); free(Bt);
}


void BuildLocalIpdgDiagHex3D(adaptive_t* adaptive, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *A) {
  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    A[n] = 0;

    // (grad phi_n, grad phi_m)_{D^e}
    for(int i=0;i<mesh->Np;++i){
      dlong base = eM*mesh->Np*mesh->Nvgeo + i;
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
      dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
      dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
      dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];    
      A[n] += JW*(dlndx*dlndx+dlndy*dlndy+dlndz*dlndz);
      A[n] += lambda*JW*B[idn]*B[idn];
    }

    for (int fM=0;fM<mesh->Nfaces;fM++) {
      // accumulate flux terms for negative and positive traces
      for(int i=0;i<mesh->Nfp;++i){
        int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

        // grab vol geofacs at surface nodes
        dlong baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
        dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
        dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
        dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
        dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
        dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
        dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];
        dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
        dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
        dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];

        // grab surface geometric factors
        dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
        dfloat nx = mesh->sgeo[base+NXID];
        dfloat ny = mesh->sgeo[base+NYID];
        dfloat nz = mesh->sgeo[base+NZID];
        dfloat wsJ = mesh->sgeo[base+WSJID];
        dfloat hinv = mesh->sgeo[base+IHID];
        
        // form negative trace terms in IPDG
        int idnM = n*mesh->Np+vidM; 
        dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
        dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
        dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
        dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
        dfloat lnM = B[idnM];
        
        dfloat penalty = adaptive->tau*hinv;     
        int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

        int bcD = 0, bcN =0;
        int bcType = 0;

        if(bc>0) bcType = adaptive->BCType[bc];          //find its type (Dirichlet/Neumann)

        // this needs to be double checked (and the code where these are used)
        if(bcType==1){ // Dirichlet
          bcD = 1;
          bcN = 0;
        } else if(bcType==2){ // Neumann
          bcD = 0;
          bcN = 1;
        }   

        A[n] += -0.5*(1+bcD)*(1-bcN)*wsJ*lnM*ndotgradlnM;  // -(ln^-, N.grad lm^-)
        A[n] += -0.5*(1+bcD)*(1-bcN)*wsJ*ndotgradlnM*lnM;  // -(N.grad ln^-, lm^-)
        A[n] += +0.5*(1+bcD)*(1-bcN)*wsJ*penalty*lnM*lnM; // +((tau/h)*ln^-,lm^-)
      }
    }
  }
}

void BuildLocalContinuousDiagHex3D(adaptive_t* adaptive, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dfloat *A) {
  for (int nz=0;nz<mesh->Nq;nz++) {
    for (int ny=0;ny<mesh->Nq;ny++) {
      for (int nx=0;nx<mesh->Nq;nx++) {
	int idn = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
	if (adaptive->mapB[idn+eM*mesh->Np]!=1) {            
	  A[idn] = 0;

	  int id = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
	  dlong base = eM*mesh->Np*mesh->Nggeo;

    
	  dfloat Grs = mesh->ggeo[base + id + G01ID*mesh->Np];
	  A[idn] += 2*Grs*mesh->D[nx+nx*mesh->Nq]*mesh->D[ny+ny*mesh->Nq];

	  dfloat Grt = mesh->ggeo[base + id + G02ID*mesh->Np];
	  A[idn] += 2*Grt*mesh->D[nx+nx*mesh->Nq]*mesh->D[nz+nz*mesh->Nq];
    
	  dfloat Gst = mesh->ggeo[base + id + G12ID*mesh->Np];
	  A[idn] += 2*Gst*mesh->D[ny+ny*mesh->Nq]*mesh->D[nz+nz*mesh->Nq];

	  for (int k=0;k<mesh->Nq;k++) {
	    int iid = k+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
	    dfloat Grr = mesh->ggeo[base + iid + G00ID*mesh->Np];
	    A[idn] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[nx+k*mesh->Nq];
	  }

	  for (int k=0;k<mesh->Nq;k++) {
	    int iid = nx+k*mesh->Nq+nz*mesh->Nq*mesh->Nq;
	    dfloat Gss = mesh->ggeo[base + iid + G11ID*mesh->Np];
	    A[idn] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[ny+k*mesh->Nq];
	  }
    
	  for (int k=0;k<mesh->Nq;k++) {
	    int iid = nx+ny*mesh->Nq+k*mesh->Nq*mesh->Nq;
	    dfloat Gtt = mesh->ggeo[base + iid + G22ID*mesh->Np];
	    A[idn] += Gtt*mesh->D[nz+k*mesh->Nq]*mesh->D[nz+k*mesh->Nq];
	  }
      
	  dfloat JW = mesh->ggeo[base + id + GWJID*mesh->Np];
	  A[idn] += JW*lambda;
	} else {
	  A[idn] = 1; //just put a 1 so A is invertable
	}
      }
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (adaptive->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (adaptive->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        A[n] += adaptive->allNeumannPenalty*adaptive->allNeumannScale*adaptive->allNeumannScale;
      } 
    }
  }
}
