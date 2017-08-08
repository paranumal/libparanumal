#include "ellipticTri2D.h"



void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the patch A matrix for element eM
void BlockJacobiPatchAx(mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A) {

  iint vbase;
  dfloat drdx;
  dfloat drdy;
  dfloat dsdx;
  dfloat dsdy;
  dfloat J;

  if (eM==-1) { //reference patch
    //equilateral triangle V = {(-1,0),(1,0),(0,sqrt(3))}
    drdx = 1.0;
    drdy = 1./sqrt(3.);
    dsdx = 0.;
    dsdy = 2./sqrt(3.);
    J = sqrt(3.)/2.;
  } else {
    vbase = eM*mesh->Nvgeo;
    drdx = mesh->vgeo[vbase+RXID];
    drdy = mesh->vgeo[vbase+RYID];
    dsdx = mesh->vgeo[vbase+SXID];
    dsdy = mesh->vgeo[vbase+SYID];
    J = mesh->vgeo[vbase+JID];
  }

  /* start with stiffness matrix  */
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){
      A[n*mesh->Np+m]  = J*lambda*mesh->MM[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*drdx*drdx*mesh->Srr[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*drdx*dsdx*mesh->Srs[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*dsdx*drdx*mesh->Ssr[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*dsdx*dsdx*mesh->Sss[n*mesh->Np+m];

      A[n*mesh->Np+m] += J*drdy*drdy*mesh->Srr[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*drdy*dsdy*mesh->Srs[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*dsdy*drdy*mesh->Ssr[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*dsdy*dsdy*mesh->Sss[n*mesh->Np+m];
    }
  }

  for (iint fM=0;fM<mesh->Nfaces;fM++) {
    iint sid;
    dfloat nx;
    dfloat ny;
    dfloat sJ;
    dfloat hinv;
    int bc=0;

    // load surface geofactors for this face
    if (eM==-1) { //reference patch
      if (fM==0) {
        nx = 0.;
        ny = -1.;
      } else if (fM==1) {
        nx = sqrt(3.)/2.;
        ny = 1./2.;
      } else if (fM==2) {
        nx = -sqrt(3.)/2.;
        ny = 1./2.;
      }
      sJ = 1.0;
      hinv = 2./sqrt(3.);
    } else {
      sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
      nx = mesh->sgeo[sid+NXID];
      ny = mesh->sgeo[sid+NYID];
      sJ = mesh->sgeo[sid+SJID];
      hinv = mesh->sgeo[sid+IHID];

      bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
    }

    dfloat penalty = tau*hinv;

    int bcD = 0, bcN =0;
    iint bcType = 0;

    if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

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
    for(iint n=0;n<mesh->Nfp;++n){
      for(iint m=0;m<mesh->Nfp;++m){
        iint nM = mesh->faceNodes[fM*mesh->Nfp+n];
        iint mM = mesh->faceNodes[fM*mesh->Nfp+m];

        // OP11 = OP11 + 0.5*( gtau*mmE )
        dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
        A[nM*mesh->Np+mM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
      }
    }

    // now add differential surface terms
    for(iint n=0;n<mesh->Nfp;++n){
      for(iint m=0;m<mesh->Np;++m){
        iint nM = mesh->faceNodes[fM*mesh->Nfp+n];

        for(iint i=0;i<mesh->Nfp;++i){
          iint iM = mesh->faceNodes[fM*mesh->Nfp+i];

          dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in

          dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+m] + dsdx*mesh->Ds[iM*mesh->Np+m];
          dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+m] + dsdy*mesh->Ds[iM*mesh->Np+m];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          A[nM*mesh->Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          A[nM*mesh->Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
        }
      }
    }

    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Nfp;++m){
        iint mM = mesh->faceNodes[fM*mesh->Nfp+m];

        for(iint i=0;i<mesh->Nfp;++i){
          iint iM = mesh->faceNodes[fM*mesh->Nfp+i];

          dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

          dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
          dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];

          // OP11 = OP11 + (- Dn1'*mmE );
          A[n*mesh->Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
          A[n*mesh->Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
        }
      }
    }
  }
}


void ellipticBuildExactBlockJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType, dfloat **patchesInvA, const char *options){


  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(iint n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Nfp;n++) {
      iint fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (iint m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;

        for (iint i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Nelements*mesh->Np*mesh->Np, sizeof(dfloat));

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){
    dfloat *patchA = patchesInvA[0] + eM*mesh->Np*mesh->Np;

    //build the patch A matrix for this element
    BlockJacobiPatchAx(mesh, basis, tau, lambda, BCType, MS, eM, patchA);

    // in place inverse (patchA points into patchesInvA[0])
    matrixInverse(mesh->Np, patchA);
  }

  free(MS);
}

void ellipticBuildApproxBlockJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   iint *Npatches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options){

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(iint n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Nfp;n++) {
      iint fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (iint m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;

        for (iint i=0;i<mesh->Np;i++){
          MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
        }

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
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

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BlockJacobiPatchAx(mesh, basis, tau, lambda, BCType, MS, -1, refPatchInvA);
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BlockJacobiPatchAx(mesh, basis, tau, lambda, BCType, MS, eM, patchA);

    iint eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    iint eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    iint eP2 = mesh->EToE[eM*mesh->Nfaces+2];

    iint fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    iint fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    iint fP2 = mesh->EToF[eM*mesh->Nfaces+2];

    if(eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interiour patch

      //hit the patch with the reference inverse
      for(iint n=0;n<mesh->Np;++n){
        for(iint m=0;m<mesh->Np;++m){
          invRefAA[n*mesh->Np+m] = 0.;
          for (iint k=0;k<mesh->Np;k++) {
            invRefAA[n*mesh->Np+m] += refPatchInvA[n*mesh->Np+k]*patchA[k*mesh->Np+m];
          }
        }
      }

      dfloat cond = matrixConditionNumber(mesh->Np,invRefAA);
      dfloat rate = (sqrt(cond)-1.)/(sqrt(cond)+1.);

      //printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < 0.7) {
        (*patchesIndex)[eM] = 0;
        refPatches++;
        continue;
      }
    }
    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA, (*Npatches)*mesh->Np*mesh->Np*sizeof(dfloat));

    matrixInverse(mesh->Np, patchA);

    //copy inverse into patchesInvA
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
        iint id = ((*Npatches)-1)*mesh->Np*mesh->Np + n*mesh->Np + m;
        (*patchesInvA)[id] = patchA[n*mesh->Np+m];
      }
    }

    (*patchesIndex)[eM] = (*Npatches)-1;
  }

  printf("saving %d full patches\n",*Npatches);
  printf("using %d reference patches\n", refPatches);

  free(patchA); free(invRefAA);
  free(MS);
}
