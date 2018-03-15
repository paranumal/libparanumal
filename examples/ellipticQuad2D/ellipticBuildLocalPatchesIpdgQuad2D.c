#include "ellipticQuad2D.h"

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the patch A matrix for element eM
void BuildLocalPatchAx(mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, iint eM, dfloat *A);


void ellipticBuildLocalPatchesIpdgQuad2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType, dfloat rateTolerance,
                                   iint *Npatches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options){

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(iint n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // build some monolithic basis arrays
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

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

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  *patchesIndex = (iint*) calloc(mesh->Nelements, sizeof(iint));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  (*Npatches) = 1;
  int refPatches = 0;


  //build a mini mesh struct for the reference patch
  mesh2D *refMesh = (mesh2D*) calloc(1,sizeof(mesh2D));
  memcpy(refMesh,mesh,sizeof(mesh2D));

  //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x =  1., V4x = -1.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y =  1.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;  refMesh->EY[0] = V1y;
  refMesh->EX[1] = V2x;  refMesh->EY[1] = V2y;
  refMesh->EX[2] = V3x;  refMesh->EY[2] = V3y;
  refMesh->EX[3] = V4x;  refMesh->EY[3] = V4y;

  refMesh->EToV = (iint*) calloc(mesh->Nverts, sizeof(iint));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;

  refMesh->EToB = (iint*) calloc(mesh->Nfaces,sizeof(iint));
  for (iint n=0;n<mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesQuad2D(refMesh, mesh->N);
  meshPhysicalNodesQuad2D(refMesh);
  meshGeometricFactorsQuad2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsQuad2D(refMesh);

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildLocalPatchAx(refMesh, basis, tau, lambda, BCType, B,Br,Bs, 0, refPatchInvA);
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildLocalPatchAx(mesh, basis, tau, lambda, BCType, B,Br,Bs, eM, patchA);

    iint eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    iint eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    iint eP2 = mesh->EToE[eM*mesh->Nfaces+2];
    iint eP3 = mesh->EToE[eM*mesh->Nfaces+3];

    iint fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    iint fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    iint fP2 = mesh->EToF[eM*mesh->Nfaces+2];
    iint fP3 = mesh->EToF[eM*mesh->Nfaces+3];

    if(eP0>=0 && eP1>=0 && eP2>=0 && eP3>=0){ //check if this is an interior patch

      refPatchInvA = *patchesInvA;

      //hit the patch with the reference inverse
      for(int n=0;n<mesh->Np;++n){
        for(int m=0;m<mesh->Np;++m){
          invRefAA[n*mesh->Np+m] = 0.;
          for (int k=0;k<mesh->Np;k++) {
            invRefAA[n*mesh->Np+m] += refPatchInvA[n*mesh->Np+k]*patchA[k*mesh->Np+m];
          }
        }
      }

      dfloat cond = matrixConditionNumber(mesh->Np,invRefAA);
      dfloat rate = (sqrt(cond)-1.)/(sqrt(cond)+1.);

      // printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = 0;
        refPatches++;
        continue;
      }
    }
    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA, (*Npatches)*mesh->Np*mesh->Np*sizeof(dfloat));

    matrixInverse(mesh->Np, patchA);

    //copy inverse into patchesInvA
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        int id = ((*Npatches)-1)*mesh->Np*mesh->Np + n*mesh->Np + m;
        (*patchesInvA)[id] = patchA[n*mesh->Np+m];
      }
    }

    (*patchesIndex)[eM] = (*Npatches)-1;
  }

  printf("saving %d full patches\n",*Npatches);
  printf("using %d reference patches\n", refPatches);

  free(refMesh);
  free(patchA); free(invRefAA);
  free(B); free(Br); free(Bs);
}


//returns the patch A matrix for element eM
void BuildLocalPatchAx(mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, iint eM, dfloat *A) {

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      A[n*mesh->Np+m] = 0;

      // (grad phi_n, grad phi_m)_{D^e}
      for(int i=0;i<mesh->Np;++i){
        iint base = eM*mesh->Np*mesh->Nvgeo + i;
        dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
        dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
        dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
        dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
        dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
        
        int idn = n*mesh->Np+i;
        int idm = m*mesh->Np+i;
        dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
        dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
        dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm];
        dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm];
        A[n*mesh->Np+m] += JW*(dlndx*dlmdx+dlndy*dlmdy);
        A[n*mesh->Np+m] += lambda*JW*B[idn]*B[idm];
      }

      for (int fM=0;fM<mesh->Nfaces;fM++) {
        // accumulate flux terms for negative and positive traces
        for(int i=0;i<mesh->Nfp;++i){
          int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

          // grab vol geofacs at surface nodes
          iint baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
          dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
          dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
          dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
          dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];

          // grab surface geometric factors
          iint base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
          dfloat nx = mesh->sgeo[base+NXID];
          dfloat ny = mesh->sgeo[base+NYID];
          dfloat wsJ = mesh->sgeo[base+WSJID];
          dfloat hinv = mesh->sgeo[base+IHID];
          
          // form negative trace terms in IPDG
          int idnM = n*mesh->Np+vidM; 
          int idmM = m*mesh->Np+vidM;

          dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
          dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
          dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
          dfloat lnM = B[idnM];

          dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];
          dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];
          dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM;
          dfloat lmM = B[idmM];
          
          dfloat penalty = tau*hinv;     
          int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

          int bcD = 0, bcN =0;
          int bcType = 0;

          if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

          // this needs to be double checked (and the code where these are used)
          if(bcType==1){ // Dirichlet
            bcD = 1;
            bcN = 0;
          } else if(bcType==2){ // Neumann
            bcD = 0;
            bcN = 1;
          }   

          A[n*mesh->Np+m] += -0.5*(1+bcD)*(1-bcN)*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
          A[n*mesh->Np+m] += -0.5*(1+bcD)*(1-bcN)*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
          A[n*mesh->Np+m] += +0.5*(1+bcD)*(1-bcN)*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)
        }
      }
    }
  }
}

