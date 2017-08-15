#include "ellipticTet3D.h"

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the patch A matrix for element eM
void BuildLocalPatchAx(mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A);

void ellipticBuildLocalPatchesIpdgTet3D(mesh3D *mesh, iint basisNp, dfloat *basis,
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

  (*Npatches) = 1;
  int refPatches = 0;

  //build a mini mesh struct for the reference patch
  mesh3D *refMesh = (mesh3D*) calloc(1,sizeof(mesh3D));
  memcpy(refMesh,mesh,sizeof(mesh3D));

  //vertices of reference patch
  dfloat V1x = -1., V2x = 1., V3x =        0., V4x = 0;
  dfloat V1y =  0., V2y = 0., V3y =  sqrt(3.), V4y = 1./sqrt(3.);
  dfloat V1z =  0., V2z = 0., V3z =        0., V4z = 2*sqrt(6.)/3.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EZ = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;  refMesh->EY[0] = V1y; refMesh->EZ[0] = V1z;
  refMesh->EX[1] = V2x;  refMesh->EY[1] = V2y; refMesh->EZ[1] = V2z;
  refMesh->EX[2] = V3x;  refMesh->EY[2] = V3y; refMesh->EZ[2] = V3z;
  refMesh->EX[3] = V4x;  refMesh->EY[3] = V4y; refMesh->EZ[3] = V4z;

  refMesh->EToV = (iint*) calloc(mesh->Nverts, sizeof(iint));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;

  refMesh->EToB = (iint*) calloc(mesh->Nfaces,sizeof(iint));
  for (iint n=0;n<mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTet3D(refMesh, mesh->N);
  meshPhysicalNodesTet3D(refMesh);
  meshGeometricFactorsTet3D(refMesh);
  meshConnectFaceNodes3D(refMesh);
  meshSurfaceGeometricFactorsTet3D(refMesh);

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  *patchesIndex = (iint*) calloc(mesh->Nelements, sizeof(iint));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildLocalPatchAx(refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  matrixInverse(mesh->Np, refPatchInvA);

  dfloat maxRate =0.;
  dfloat maxCond =0.;

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildLocalPatchAx(mesh, basis, tau, lambda, BCType, MS, eM, patchA);

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
      maxRate = mymax(rate,maxRate);
      maxCond = mymax(cond,maxCond);

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
  printf("Max condition number = %g, and slowest CG convergence rate = %g\n", maxCond, maxRate);


  free(refMesh);
  free(patchA); free(invRefAA);
  free(MS);
}

//returns the patch A matrix for element eM
void BuildLocalPatchAx(mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A) {

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
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){

      A[n*mesh->Np+m]  = J*lambda*mesh->MM[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G00*mesh->Srr[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G01*mesh->Srs[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G02*mesh->Srt[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G10*mesh->Ssr[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G11*mesh->Sss[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G12*mesh->Sst[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G20*mesh->Str[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G21*mesh->Sts[n*mesh->Np+m];
      A[n*mesh->Np+m] += J*G22*mesh->Stt[n*mesh->Np+m];
    }
  }

  for (iint fM=0;fM<mesh->Nfaces;fM++) {
    // load surface geofactors for this face
    iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat nz = mesh->sgeo[sid+NZID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat hinv = mesh->sgeo[sid+IHID];

    int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

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

          dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+m] + dsdx*mesh->Ds[iM*mesh->Np+m] + dtdx*mesh->Dt[iM*mesh->Np+m];
          dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+m] + dsdy*mesh->Ds[iM*mesh->Np+m] + dtdy*mesh->Dt[iM*mesh->Np+m];
          dfloat DzMim = drdz*mesh->Dr[iM*mesh->Np+m] + dsdz*mesh->Ds[iM*mesh->Np+m] + dtdz*mesh->Dt[iM*mesh->Np+m];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          A[nM*mesh->Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          A[nM*mesh->Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
          A[nM*mesh->Np+m] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;
        }
      }
    }

    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Nfp;++m){
        iint mM = mesh->faceNodes[fM*mesh->Nfp+m];

        for(iint i=0;i<mesh->Nfp;++i){
          iint iM = mesh->faceNodes[fM*mesh->Nfp+i];

          dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

          dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n] + dtdx*mesh->Dt[iM*mesh->Np+n];
          dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n] + dtdy*mesh->Dt[iM*mesh->Np+n];
          dfloat DzMin = drdz*mesh->Dr[iM*mesh->Np+n] + dsdz*mesh->Ds[iM*mesh->Np+n] + dtdz*mesh->Dt[iM*mesh->Np+n];

          // OP11 = OP11 + (- Dn1'*mmE );
          A[n*mesh->Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
          A[n*mesh->Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
          A[n*mesh->Np+mM] +=  -0.5*nz*(1+bcD)*(1-bcN)*DzMin*MSfim;
        }
      }
    }
  }
}

