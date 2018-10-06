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

#include "elliptic.h"

void BuildLocalIpdgBBDiagTri2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dlong eM, dfloat *A);
void BuildLocalIpdgDiagTri2D (elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dlong eM, dfloat *A);
void BuildLocalIpdgDiagQuad2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dlong eM, dfloat *A);
void BuildLocalIpdgDiagQuad3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dlong eM, dfloat *A);
void BuildLocalIpdgDiagTet3D (elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dlong eM, dfloat *A);
void BuildLocalIpdgDiagHex3D (elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *A);

void BuildLocalContinuousDiagTri2D (elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *A);
void BuildLocalContinuousDiagQuad2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *A);
void BuildLocalContinuousDiagQuad3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *A);
void BuildLocalContinuousDiagTet3D (elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *A);
void BuildLocalContinuousDiagHex3D (elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dfloat *A);

void ellipticBuildJacobi(elliptic_t* elliptic, dfloat lambda, dfloat **invDiagA){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

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

  if (elliptic->elementType==QUADRILATERALS) {
    int mode = 0;
    for(int nj=0;nj<mesh->N+1;++nj){
      for(int ni=0;ni<mesh->N+1;++ni){

        int node = 0;

        for(int j=0;j<mesh->N+1;++j){
          for(int i=0;i<mesh->N+1;++i){

            if(nj==j && ni==i)
              B[mode*mesh->Np+node] = 1;
            if(nj==j)
              Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
            if(ni==i)
              Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
            
            ++node;
          }
        }
        ++mode;
      }
    }
  }

  if (elliptic->elementType==HEXAHEDRA) {
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
    switch(elliptic->elementType){
      case TRIANGLES: 
        if (options.compareArgs("BASIS","BERN")) {
          #pragma omp parallel for 
          for(dlong eM=0;eM<mesh->Nelements;++eM)
            BuildLocalIpdgBBDiagTri2D(elliptic, mesh, lambda, MS, eM, diagA + eM*mesh->Np);   
        } else {
          #pragma omp parallel for 
          for(dlong eM=0;eM<mesh->Nelements;++eM)
            BuildLocalIpdgDiagTri2D(elliptic, mesh, lambda, MS, eM, diagA + eM*mesh->Np); 
        } 
        break;
      case QUADRILATERALS:
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM)
          BuildLocalIpdgDiagQuad2D(elliptic, mesh, lambda, MS, B, Br, Bs, eM, diagA + eM*mesh->Np);
        break;
      case TETRAHEDRA:
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM)
          BuildLocalIpdgDiagTet3D(elliptic, mesh, lambda, MS, eM, diagA + eM*mesh->Np); 
        break;
      case HEXAHEDRA:
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM)
          BuildLocalIpdgDiagHex3D(elliptic, mesh, lambda, MS, B, Br, Bs, Bt, eM, diagA + eM*mesh->Np); 
        break;
    }
  } else if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    switch(elliptic->elementType){
      case TRIANGLES: 
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM)
          BuildLocalContinuousDiagTri2D(elliptic, mesh, lambda, eM, diagA + eM*mesh->Np);
        break;
      case QUADRILATERALS:{
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM){
          if(elliptic->dim==2)
            BuildLocalContinuousDiagQuad2D(elliptic, mesh, lambda, eM, B, Br, Bs, diagA + eM*mesh->Np);
          if(elliptic->dim==3)
            BuildLocalContinuousDiagQuad3D(elliptic, mesh, lambda, eM, B, Br, Bs, diagA + eM*mesh->Np);
          }
        }break;
      case TETRAHEDRA:
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM)
          BuildLocalContinuousDiagTet3D(elliptic, mesh, lambda, eM, diagA + eM*mesh->Np);
        break;
      case HEXAHEDRA:
        #pragma omp parallel for 
        for(dlong eM=0;eM<mesh->Nelements;++eM)
          BuildLocalContinuousDiagHex3D(elliptic, mesh, lambda, eM, B, Br, Bs, Bt, diagA + eM*mesh->Np);
        break;
    }
  }

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) 
    ogsGatherScatter(diagA, ogsDfloat, ogsAdd, elliptic->ogs);
    
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(mesh->rank==0) printf("done.\n");

  free(diagA);
  free(MS);
  free(B); free(Br); free(Bs); free(Bt);
}

void BuildLocalIpdgDiagTri2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dlong eM, dfloat *A) {

  dlong vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];
  dfloat J = mesh->vgeo[vbase+JID];

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    A[n]  = J*lambda*mesh->MM[n*mesh->Np+n];
    A[n] += J*drdx*drdx*mesh->Srr[n*mesh->Np+n];
    A[n] += J*drdx*dsdx*mesh->Srs[n*mesh->Np+n];
    A[n] += J*dsdx*drdx*mesh->Ssr[n*mesh->Np+n];
    A[n] += J*dsdx*dsdx*mesh->Sss[n*mesh->Np+n];
    A[n] += J*drdy*drdy*mesh->Srr[n*mesh->Np+n];
    A[n] += J*drdy*dsdy*mesh->Srs[n*mesh->Np+n];
    A[n] += J*dsdy*drdy*mesh->Ssr[n*mesh->Np+n];
    A[n] += J*dsdy*dsdy*mesh->Sss[n*mesh->Np+n];
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      A[n] += elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
    }
  }

  for (int fM=0;fM<mesh->Nfaces;fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat hinv = mesh->sgeo[sid+IHID];

    int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

    dfloat penalty = elliptic->tau*hinv;

    int bcD = 0, bcN =0;
    int bcType = 0;

    if(bc>0) bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)

    // this needs to be double checked (and the code where these are used)
    if(bcType==1){ // Dirichlet
      bcD = 1;
      bcN = 0;
    } else if(bcType==2){ // Neumann
      bcD = 0;
      bcN = 1;
    }

    // mass matrix for this face
    dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

    // penalty term just involves face nodes
    for(int n=0;n<mesh->Nfp;++n){
      int nM = mesh->faceNodes[fM*mesh->Nfp+n];
      
      for(int m=0;m<mesh->Nfp;++m){
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];
        if (mM == nM) {
          // OP11 = OP11 + 0.5*( gtau*mmE )
          dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
          A[nM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
        }
      }
    }

    // now add differential surface terms
    for(int n=0;n<mesh->Nfp;++n){
      int nM = mesh->faceNodes[fM*mesh->Nfp+n];

      for(int i=0;i<mesh->Nfp;++i){
        int iM = mesh->faceNodes[fM*mesh->Nfp+i];

        dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in

        dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+nM] + dsdx*mesh->Ds[iM*mesh->Np+nM];
        dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+nM] + dsdy*mesh->Ds[iM*mesh->Np+nM];

        // OP11 = OP11 + 0.5*( - mmE*Dn1)
        A[nM] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
        A[nM] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
      }
    }

    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Nfp;++m){
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];

        if (mM==n) {
          for(int i=0;i<mesh->Nfp;++i){
            int iM = mesh->faceNodes[fM*mesh->Nfp+i];

            dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];

            // OP11 = OP11 + (- Dn1'*mmE );
            A[n] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            A[n] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }
    }
  }
}

void BuildLocalIpdgPatchAxTri2D(elliptic_t* elliptic, mesh_t* mesh, int basisNp, dfloat *basis, dfloat lambda,
                        dfloat *MS, dlong eM, dfloat *A);

//generate the BB diagonal by extracting it from the transformed patch
void BuildLocalIpdgBBDiagTri2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dlong eM, dfloat *A) {

  dfloat *patchA = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  int basisNp = mesh->Np;
  dfloat *basis = mesh->VB;

  BuildLocalIpdgPatchAxTri2D(elliptic, mesh, basisNp, basis, lambda, MS, eM, patchA);

  for(int n=0;n<mesh->Np;++n) {
    A[n] = patchA[n*mesh->Np+n]; //store the diagonal entry
  }
  free(patchA);
}

//returns the continuous C0 patch A matrix for element eM
void BuildLocalContinuousDiagTri2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *A) {

  dlong gbase = eM*mesh->Nggeo;
  dfloat Grr = mesh->ggeo[gbase + G00ID];
  dfloat Grs = mesh->ggeo[gbase + G01ID];
  dfloat Gss = mesh->ggeo[gbase + G11ID];
  dfloat J   = mesh->ggeo[gbase + GWJID];

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    if (elliptic->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
      A[n] = J*lambda*mesh->MM[n+n*mesh->Np];
      A[n] += Grr*mesh->Srr[n+n*mesh->Np];
      A[n] += Grs*mesh->Srs[n+n*mesh->Np];
      A[n] += Grs*mesh->Ssr[n+n*mesh->Np];
      A[n] += Gss*mesh->Sss[n+n*mesh->Np];
    } else {
      A[n] = 1; //just put a 1 so A is invertable
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (elliptic->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        A[n] += elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
      } 
    }
  }
}

void BuildLocalIpdgDiagQuad2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dlong eM, dfloat *A) {
  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    A[n] = 0;

    // (grad phi_n, grad phi_m)_{D^e}
    for(int i=0;i<mesh->Np;++i){
      dlong base = eM*mesh->Np*mesh->Nvgeo + i;
      dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
      dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
      dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
      dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
      dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
      
      int idn = n*mesh->Np+i;
      dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
      dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
      A[n] += JW*(dlndx*dlndx+dlndy*dlndy);
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
        dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
        dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];

        // grab surface geometric factors
        dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
        dfloat nx = mesh->sgeo[base+NXID];
        dfloat ny = mesh->sgeo[base+NYID];
        dfloat wsJ = mesh->sgeo[base+WSJID];
        dfloat hinv = mesh->sgeo[base+IHID];
        
        // form negative trace terms in IPDG
        int idnM = n*mesh->Np+vidM; 

        dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
        dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
        dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
        dfloat lnM = B[idnM];
        
        dfloat penalty = elliptic->tau*hinv;     
        int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

        int bcD = 0, bcN =0;
        int bcType = 0;

        if(bc>0) bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)

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

void BuildLocalContinuousDiagQuad2D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, 
                                    dfloat *B, dfloat *Br, dfloat* Bs, dfloat *A) {

  for (int ny=0;ny<mesh->Nq;ny++) {
    for (int nx=0;nx<mesh->Nq;nx++) {
      int iid = nx+ny*mesh->Nq;
      if (elliptic->mapB[nx+ny*mesh->Nq+eM*mesh->Np]!=1) {    
        A[iid] = 0;

        for (int k=0;k<mesh->Nq;k++) {
          int id = k+ny*mesh->Nq;
          dfloat Grr = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];
          A[iid] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[nx+k*mesh->Nq];
        }

        for (int k=0;k<mesh->Nq;k++) {
          int id = nx+k*mesh->Nq;
          dfloat Gss = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];
          A[iid] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[ny+k*mesh->Nq];
        }

        int id = nx+ny*mesh->Nq;
        dfloat Grs = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
        A[iid] += 2*Grs*mesh->D[nx+nx*mesh->Nq]*mesh->D[ny+ny*mesh->Nq];

        dfloat JW = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + GWJID*mesh->Np];
        A[iid] += JW*lambda;

      } else {
        A[iid] = 1; //just put a 1 so A is invertable
      }
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (elliptic->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        A[n] += elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
      } 
    }
  }
}


void BuildLocalIpdgDiagQuad3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dlong eM, dfloat *A) {
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
      dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx;
      dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy;
      dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz;
      A[n] += JW*(dlndx*dlndx + dlndy*dlndy + dlndz*dlndz);
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

        dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM;
        dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM;
        dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM;
        dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM ;
        dfloat lnM = B[idnM];
        
        dfloat penalty = elliptic->tau*hinv;     
       
        A[n] += -0.5*wsJ*lnM*ndotgradlnM;  // -(ln^-, N.grad lm^-)
        A[n] += -0.5*wsJ*ndotgradlnM*lnM;  // -(N.grad ln^-, lm^-)
        A[n] += +0.5*wsJ*penalty*lnM*lnM; // +((tau/h)*ln^-,lm^-)
      }
    }
  }
}

void BuildLocalContinuousDiagQuad3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, 
                                    dfloat *B, dfloat *Br, dfloat* Bs, dfloat *A) {

  for (int ny=0;ny<mesh->Nq;ny++) {
    for (int nx=0;nx<mesh->Nq;nx++) {
      int iid = nx+ny*mesh->Nq;
      A[iid] = 0;

      for (int k=0;k<mesh->Nq;k++) {
        int id = k+ny*mesh->Nq;
        dfloat Grr = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];
        A[iid] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[nx+k*mesh->Nq];
      }

      for (int k=0;k<mesh->Nq;k++) {
        int id = nx+k*mesh->Nq;
        dfloat Gss = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];
        A[iid] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[ny+k*mesh->Nq];
      }

      int id = nx+ny*mesh->Nq;
      dfloat Grs = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
      A[iid] += 2*Grs*mesh->D[nx+nx*mesh->Nq]*mesh->D[ny+ny*mesh->Nq];

      // id = nx+ny*mesh->Nq;
      // dfloat Grt = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
      // A[iid] += 2*Grt*mesh->D[nx+nx*mesh->Nq];

      // id = nx+ny*mesh->Nq;
      // dfloat Gst = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
      // A[iid] += 2*Gst*mesh->D[ny+ny*mesh->Nq];
      
      // dfloat Gtt = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G22ID*mesh->Np];
      // A[iid] += Gtt;
      


      dfloat JW  = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + GWJID*mesh->Np];
      A[iid] += JW*lambda;
    }
  }
}



void BuildLocalIpdgDiagTet3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dlong eM, dfloat *A) {
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
  dfloat J = mesh->vgeo[vbase+JID];

  dfloat G00 = drdx*drdx + drdy*drdy + drdz*drdz;
  dfloat G01 = drdx*dsdx + drdy*dsdy + drdz*dsdz;
  dfloat G02 = drdx*dtdx + drdy*dtdy + drdz*dtdz;

  dfloat G10 = dsdx*drdx + dsdy*drdy + dsdz*drdz;
  dfloat G11 = dsdx*dsdx + dsdy*dsdy + dsdz*dsdz;
  dfloat G12 = dsdx*dtdx + dsdy*dtdy + dsdz*dtdz;

  dfloat G20 = dtdx*drdx + dtdy*drdy + dtdz*drdz;
  dfloat G21 = dtdx*dsdx + dtdy*dsdy + dtdz*dsdz;
  dfloat G22 = dtdx*dtdx + dtdy*dtdy + dtdz*dtdz;


  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    A[n]  = J*lambda*mesh->MM[n*mesh->Np+n];
    A[n] += J*G00*mesh->Srr[n*mesh->Np+n];
    A[n] += J*G01*mesh->Srs[n*mesh->Np+n];
    A[n] += J*G02*mesh->Srt[n*mesh->Np+n];
    A[n] += J*G10*mesh->Ssr[n*mesh->Np+n];
    A[n] += J*G11*mesh->Sss[n*mesh->Np+n];
    A[n] += J*G12*mesh->Sst[n*mesh->Np+n];
    A[n] += J*G20*mesh->Str[n*mesh->Np+n];
    A[n] += J*G21*mesh->Sts[n*mesh->Np+n];
    A[n] += J*G22*mesh->Stt[n*mesh->Np+n];
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      A[n] += elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
    }
  }

  for (int fM=0;fM<mesh->Nfaces;fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat nz = mesh->sgeo[sid+NZID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat hinv = mesh->sgeo[sid+IHID];

    int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

    dfloat penalty = elliptic->tau*hinv;

    int bcD = 0, bcN =0;
    int bcType = 0;

    if(bc>0) bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)

    // this needs to be double checked (and the code where these are used)
    if(bcType==1){ // Dirichlet
      bcD = 1;
      bcN = 0;
    } else if(bcType==2){ // Neumann
      bcD = 0;
      bcN = 1;
    }

    // mass matrix for this face
    dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

    // penalty term just involves face nodes
    for(int n=0;n<mesh->Nfp;++n){
      for(int m=0;m<mesh->Nfp;++m){
        int nM = mesh->faceNodes[fM*mesh->Nfp+n];
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];

        if (mM==nM) {
          // OP11 = OP11 + 0.5*( gtau*mmE )
          dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
          A[nM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
        }
      }
    }

    // now add differential surface terms
    for(int n=0;n<mesh->Nfp;++n){
      int nM = mesh->faceNodes[fM*mesh->Nfp+n];

      for(int i=0;i<mesh->Nfp;++i){
        int iM = mesh->faceNodes[fM*mesh->Nfp+i];

        dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in

        dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+nM] + dsdx*mesh->Ds[iM*mesh->Np+nM] + dtdx*mesh->Dt[iM*mesh->Np+nM];
        dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+nM] + dsdy*mesh->Ds[iM*mesh->Np+nM] + dtdy*mesh->Dt[iM*mesh->Np+nM];
        dfloat DzMim = drdz*mesh->Dr[iM*mesh->Np+nM] + dsdz*mesh->Ds[iM*mesh->Np+nM] + dtdz*mesh->Dt[iM*mesh->Np+nM];

        // OP11 = OP11 + 0.5*( - mmE*Dn1)
        A[nM] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
        A[nM] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
        A[nM] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;
      }
    }

    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Nfp;++m){
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];

        if (mM==n) {
          for(int i=0;i<mesh->Nfp;++i){
            int iM = mesh->faceNodes[fM*mesh->Nfp+i];

            dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n] + dtdx*mesh->Dt[iM*mesh->Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n] + dtdy*mesh->Dt[iM*mesh->Np+n];
            dfloat DzMin = drdz*mesh->Dr[iM*mesh->Np+n] + dsdz*mesh->Ds[iM*mesh->Np+n] + dtdz*mesh->Dt[iM*mesh->Np+n];

            // OP11 = OP11 + (- Dn1'*mmE );
            A[n] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            A[n] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
            A[n] +=  -0.5*nz*(1+bcD)*(1-bcN)*DzMin*MSfim;
          }
        }
      }
    }
  }
}

void BuildLocalContinuousDiagTet3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *A) {
  dlong gbase = eM*mesh->Nggeo;
  dfloat Grr = mesh->ggeo[gbase + G00ID];
  dfloat Grs = mesh->ggeo[gbase + G01ID];
  dfloat Grt = mesh->ggeo[gbase + G02ID];
  dfloat Gss = mesh->ggeo[gbase + G11ID];
  dfloat Gst = mesh->ggeo[gbase + G12ID];
  dfloat Gtt = mesh->ggeo[gbase + G22ID];
  dfloat J   = mesh->ggeo[gbase + GWJID];

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    if (elliptic->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
      A[n] = J*lambda*mesh->MM[n+n*mesh->Np];
      A[n] += Grr*mesh->Srr[n+n*mesh->Np];
      A[n] += Grs*mesh->Srs[n+n*mesh->Np];
      A[n] += Grt*mesh->Srt[n+n*mesh->Np];
      A[n] += Grs*mesh->Ssr[n+n*mesh->Np];
      A[n] += Gss*mesh->Sss[n+n*mesh->Np];
      A[n] += Gst*mesh->Sst[n+n*mesh->Np];
      A[n] += Grt*mesh->Str[n+n*mesh->Np];
      A[n] += Gst*mesh->Sts[n+n*mesh->Np];
      A[n] += Gtt*mesh->Stt[n+n*mesh->Np];
    } else {
      A[n] = 1; //just put a 1 so A is invertable
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (elliptic->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        A[n] += elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
      } 
    }
  }
}

void BuildLocalIpdgDiagHex3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dfloat *MS, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *A) {
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
        
        dfloat penalty = elliptic->tau*hinv;     
        int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

        int bcD = 0, bcN =0;
        int bcType = 0;

        if(bc>0) bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)

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

void BuildLocalContinuousDiagHex3D(elliptic_t* elliptic, mesh_t *mesh, dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dfloat *A) {
  for (int nz=0;nz<mesh->Nq;nz++) {
  for (int ny=0;ny<mesh->Nq;ny++) {
  for (int nx=0;nx<mesh->Nq;nx++) {
    int idn = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
    if (elliptic->mapB[idn+eM*mesh->Np]!=1) {            
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
  if (elliptic->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (elliptic->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        A[n] += elliptic->allNeumannPenalty*elliptic->allNeumannScale*elliptic->allNeumannScale;
      } 
    }
  }
}
