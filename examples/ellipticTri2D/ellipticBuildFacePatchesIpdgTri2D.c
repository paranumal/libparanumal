#include "ellipticTri2D.h"

void matrixInverse(int N, dfloat *A);

dfloat matrixConditionNumber(int N, dfloat *A);

void BuildFacePatchAx(mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint face, dfloat *A);


void ellipticBuildFacePatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
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

  //We need the halo element's EToB flags to make the patch matrices
  if (mesh->totalHaloPairs) {
    mesh->EToB = (int *) realloc(mesh->EToB,(mesh->Nelements+mesh->totalHaloPairs)*sizeof(int));
    iint *idSendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, sizeof(int), mesh->EToB, idSendBuffer, mesh->EToB + mesh->Nelements);
    free(idSendBuffer);
  }

  int NpatchElements = 2;
  int patchNp = mesh->Np*NpatchElements;

  //build a list of all face pairs
  mesh->NfacePairs=0;
  for (iint eM=0; eM<mesh->Nelements;eM++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      iint eP = mesh->EToE[eM*mesh->Nfaces+f];

      if (eM<eP) mesh->NfacePairs++;
    }
  }

  mesh->EToFPairs = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfaces,sizeof(iint));
  mesh->FPairsToE = (iint *) calloc(2*mesh->NfacePairs,sizeof(iint));
  mesh->FPairsToF = (int *) calloc(2*mesh->NfacePairs,sizeof(int));

  //fill with -1
  for (iint n=0;n<(mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfaces;n++) {
    mesh->EToFPairs[n] = -1;
  }

  iint cnt=0;
  for (iint eM=0; eM<mesh->Nelements;eM++) {
    for (int fM=0;fM<mesh->Nfaces;fM++) {
      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];

      if (eM<eP) {
        mesh->FPairsToE[2*cnt+0] = eM;
        mesh->FPairsToE[2*cnt+1] = eP;

        int fP = mesh->EToF[eM*mesh->Nfaces+fM];

        mesh->FPairsToF[2*cnt+0] = fM;
        mesh->FPairsToF[2*cnt+1] = fP;

        mesh->EToFPairs[mesh->Nfaces*eM+fM] = 2*cnt+0;
        mesh->EToFPairs[mesh->Nfaces*eP+fP] = 2*cnt+1;

        cnt++;
      }
    }
  }

  //build a mini mesh struct for the reference patch
  mesh2D *refMesh = (mesh2D*) calloc(1,sizeof(mesh2D));
  memcpy(refMesh,mesh,sizeof(mesh2D));

  //vertices of reference patch
  dfloat V1x = -1., V2x = 1., V3x =        0., V4x =        0.;
  dfloat V1y =  0., V2y = 0., V3y =  sqrt(3.), V4y = -sqrt(3.);

  refMesh->Nelements = NpatchElements;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));

  refMesh->EX[0*mesh->Nverts+0] = V1x;  refMesh->EY[0*mesh->Nverts+0] = V1y;
  refMesh->EX[0*mesh->Nverts+1] = V2x;  refMesh->EY[0*mesh->Nverts+1] = V2y;
  refMesh->EX[0*mesh->Nverts+2] = V3x;  refMesh->EY[0*mesh->Nverts+2] = V3y;

  refMesh->EX[1*mesh->Nverts+0] = V2x;  refMesh->EY[1*mesh->Nverts+0] = V2y;
  refMesh->EX[1*mesh->Nverts+1] = V1x;  refMesh->EY[1*mesh->Nverts+1] = V1y;
  refMesh->EX[1*mesh->Nverts+2] = V4x;  refMesh->EY[1*mesh->Nverts+2] = V4y;

  refMesh->EToV = (iint*) calloc(NpatchElements*mesh->Nverts, sizeof(iint));

  refMesh->EToV[0*mesh->Nverts+0] = 0;
  refMesh->EToV[0*mesh->Nverts+1] = 1;
  refMesh->EToV[0*mesh->Nverts+2] = 2;

  refMesh->EToV[1*mesh->Nverts+0] = 1;
  refMesh->EToV[1*mesh->Nverts+1] = 0;
  refMesh->EToV[1*mesh->Nverts+2] = 3;

  refMesh->EToB = (iint*) calloc(NpatchElements*mesh->Nfaces,sizeof(iint));
  for (iint n=0;n<NpatchElements*mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTri2D(refMesh, mesh->N);
  meshPhysicalNodesTri2D(refMesh);
  meshGeometricFactorsTri2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsTri2D(refMesh);

  //build a list of all face pairs
  refMesh->NfacePairs=1;

  refMesh->EToFPairs = (iint *) calloc(2*mesh->Nfaces,sizeof(iint));
  refMesh->FPairsToE = (iint *) calloc(2,sizeof(iint));
  refMesh->FPairsToF = (int *)  calloc(2,sizeof(int));

  //fill with -1
  for (iint n=0;n<2*mesh->Nfaces;n++)  refMesh->EToFPairs[n] = -1;

  refMesh->FPairsToE[0] = 0;
  refMesh->FPairsToE[1] = 1;
  refMesh->FPairsToF[0] = 0;
  refMesh->FPairsToF[1] = 0;
  refMesh->EToFPairs[0] = 0;
  refMesh->EToFPairs[refMesh->Nfaces] = 1;



  iint Nperm = pow(mesh->Nfaces,2);//all possible configuration of neighbours
  (*Npatches) = Nperm;
  int refPatches = 0;

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(Nperm*patchNp*patchNp, sizeof(dfloat));
  *patchesIndex = (iint*) calloc(mesh->NfacePairs, sizeof(iint));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));


  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildFacePatchAx(refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  matrixInverse(patchNp, refPatchInvA);

  //store the permutations of the reference patch
  int blkCounter = 1;
  int *permIndex = (int*) calloc(patchNp, sizeof(int));
  for(iint blk=1;blk<Nperm;++blk){
    iint f0 = blk%mesh->Nfaces;
    iint f1 = (blk/mesh->Nfaces)%mesh->Nfaces;

    iint r0 = (f0==0) ? 0: mesh->Nfaces-f0;
    iint r1 = (f1==0) ? 0: mesh->Nfaces-f1;

    for(iint n=0;n<mesh->Np;++n){
      permIndex[n+0*mesh->Np] = 0*mesh->Np + mesh->rmapP[r0*mesh->Np+n];
      permIndex[n+1*mesh->Np] = 1*mesh->Np + mesh->rmapP[r1*mesh->Np+n];
    }

    for(iint n=0;n<patchNp;++n){
      for(iint m=0;m<patchNp;++m){
        iint pn = permIndex[n];
        iint pm = permIndex[m];
        (*patchesInvA)[blkCounter*patchNp*patchNp + n*patchNp + m] = refPatchInvA[pn*patchNp+pm]; // maybe need to switch map
      }
    }
    ++blkCounter;
  }

  // loop over all elements
  for(iint face=0;face<mesh->NfacePairs;++face){

    //build the patch A matrix for this element
    BuildFacePatchAx(mesh, basis, tau, lambda, BCType, MS, face, patchA);

    iint eM = mesh->FPairsToE[2*face+0];
    iint eP = mesh->FPairsToE[2*face+1];
    iint fM = mesh->FPairsToF[2*face+0];
    iint fP = mesh->FPairsToF[2*face+1];

    if (eP >=0) {
      iint eM0 = mesh->EToE[eM*mesh->Nfaces+0];
      iint eM1 = mesh->EToE[eM*mesh->Nfaces+1];
      iint eM2 = mesh->EToE[eM*mesh->Nfaces+2];

      iint eP0 = mesh->EToE[eP*mesh->Nfaces+0];
      iint eP1 = mesh->EToE[eP*mesh->Nfaces+1];
      iint eP2 = mesh->EToE[eP*mesh->Nfaces+2];

      if(eM0>=0 && eM1>=0 && eM2>=0 &&
          eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interiour patch
        iint blk = fM + mesh->Nfaces*fP;

        refPatchInvA = *patchesInvA + blk*patchNp*patchNp;

        //hit the patch with the reference inverse
        for(iint n=0;n<patchNp;++n){
          for(iint m=0;m<patchNp;++m){
            invRefAA[n*patchNp+m] = 0.;
            for (iint k=0;k<patchNp;k++) {
              invRefAA[n*patchNp+m] += refPatchInvA[n*patchNp+k]*patchA[k*patchNp+m];
            }
          }
        }

        dfloat cond = matrixConditionNumber(patchNp,invRefAA);
        dfloat rate = (sqrt(cond)-1.)/(sqrt(cond)+1.);

        //printf("Face %d's conditioned patch reports cond = %g and rate = %g \n", face, cond, rate);

        if (rate < rateTolerance) {
          (*patchesIndex)[face] = blk;
          refPatches++;
          continue;
        }
      }
    }

    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA, (*Npatches)*patchNp*patchNp*sizeof(dfloat));

    matrixInverse(patchNp, patchA);

    //copy inverse into patchesInvA
    for(iint n=0;n<patchNp;++n){
      for(iint m=0;m<patchNp;++m){
        iint id = ((*Npatches)-1)*patchNp*patchNp + n*patchNp + m;
        (*patchesInvA)[id] = patchA[n*patchNp+m];
      }
    }

    (*patchesIndex)[face] = (*Npatches)-1;
  }

  printf("saving %d full patches\n",*Npatches);
  printf("using %d reference patches\n", refPatches);

  free(patchA); free(invRefAA);
  free(MS);
}

void BuildFacePatchAx(mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint face, dfloat *A) {

  int NpatchElements = 2;
  int patchNp = NpatchElements*mesh->Np;

  // Extract patches
  // B  a
  // a' C

  //zero out the matrix
  for (int n=0;n<patchNp*patchNp;n++) A[n] = 0.;

  //make sure the diagonal is at least identity
  for (int n=0;n<patchNp;n++) A[n*patchNp+n] = 1.;

  //start with diagonals
  for(iint N=0;N<NpatchElements;++N){
    //element number
    iint e = mesh->FPairsToE[2*face+N];

    iint vbase = e*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
        int id = N*mesh->Np*patchNp + n*patchNp + N*mesh->Np + m;

        A[id]  = J*lambda*mesh->MM[n*mesh->Np+m];
        A[id] += J*drdx*drdx*mesh->Srr[n*mesh->Np+m];
        A[id] += J*drdx*dsdx*mesh->Srs[n*mesh->Np+m];
        A[id] += J*dsdx*drdx*mesh->Ssr[n*mesh->Np+m];
        A[id] += J*dsdx*dsdx*mesh->Sss[n*mesh->Np+m];

        A[id] += J*drdy*drdy*mesh->Srr[n*mesh->Np+m];
        A[id] += J*drdy*dsdy*mesh->Srs[n*mesh->Np+m];
        A[id] += J*dsdy*drdy*mesh->Ssr[n*mesh->Np+m];
        A[id] += J*dsdy*dsdy*mesh->Sss[n*mesh->Np+m];
      }
    }

    for (iint fM=0;fM<mesh->Nfaces;fM++) {
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(e*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];

      int bc = mesh->EToB[fM+mesh->Nfaces*e]; //raw boundary flag

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
          int nM = mesh->faceNodes[fM*mesh->Nfp+n];
          int mM = mesh->faceNodes[fM*mesh->Nfp+m];
          int id = N*mesh->Np*patchNp + nM*patchNp + N*mesh->Np + mM;

          // OP11 = OP11 + 0.5*( gtau*mmE )
          dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
          A[id] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
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

            int id = N*mesh->Np*patchNp + nM*patchNp + N*mesh->Np + m;

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            A[id] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            A[id] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
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

            int id = N*mesh->Np*patchNp + n*patchNp + N*mesh->Np + mM;

            // OP11 = OP11 + (- Dn1'*mmE );
            A[id] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            A[id] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }
    }
  }

  //now the off-diagonal
  iint eM = mesh->FPairsToE[2*face+0];
  iint eP = mesh->FPairsToE[2*face+1];
  int fM = mesh->FPairsToF[2*face+0];

  iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
  dfloat nx = mesh->sgeo[sid+NXID];
  dfloat ny = mesh->sgeo[sid+NYID];
  dfloat sJ = mesh->sgeo[sid+SJID];
  dfloat hinv = mesh->sgeo[sid+IHID];

  iint vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];
  dfloat J = mesh->vgeo[vbase+JID];

  iint vbaseP = eP*mesh->Nvgeo;
  dfloat drdxP = mesh->vgeo[vbaseP+RXID];
  dfloat drdyP = mesh->vgeo[vbaseP+RYID];
  dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
  dfloat dsdyP = mesh->vgeo[vbaseP+SYID];

  dfloat penalty = tau*hinv;

  // mass matrix for this face
  dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

  // penalty term just involves face nodes
  for(iint n=0;n<mesh->Nfp;++n){
    for(iint m=0;m<mesh->Nfp;++m){
      iint nM = mesh->faceNodes[fM*mesh->Nfp+n];
      iint mM = mesh->faceNodes[fM*mesh->Nfp+m];

      dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];

      // neighbor penalty term
      iint idM = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+m;
      int mP = mesh->vmapP[idM]%mesh->Np;

      int id = nM*patchNp + mesh->Np + mP;

      // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
      A[id] += -0.5*penalty*MSfnm;
    }
  }

  // now add differential surface terms
  for(iint n=0;n<mesh->Nfp;++n){
    for(iint m=0;m<mesh->Np;++m){
      int nM = mesh->faceNodes[fM*mesh->Nfp+n];

      for(iint i=0;i<mesh->Nfp;++i){
        int iM = mesh->faceNodes[fM*mesh->Nfp+i];
        int iP = mesh->vmapP[i + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;

        dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in

        dfloat DxPim = drdxP*mesh->Dr[iP*mesh->Np+m] + dsdxP*mesh->Ds[iP*mesh->Np+m];
        dfloat DyPim = drdyP*mesh->Dr[iP*mesh->Np+m] + dsdyP*mesh->Ds[iP*mesh->Np+m];

        int id = nM*patchNp + mesh->Np + m;

        //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
        A[id] += -0.5*nx*MSfni*DxPim;
        A[id] += -0.5*ny*MSfni*DyPim;
      }
    }
  }

  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Nfp;++m){
      int mM = mesh->faceNodes[fM*mesh->Nfp+m];
      int mP = mesh->vmapP[m + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;

      for(iint i=0;i<mesh->Nfp;++i){
        iint iM = mesh->faceNodes[fM*mesh->Nfp+i];

        dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

        dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
        dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];

        int id = n*patchNp + mesh->Np + mP;

        //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
        A[id] +=  +0.5*nx*DxMin*MSfim;
        A[id] +=  +0.5*ny*DyMin*MSfim;
      }
    }
  }

  //write the transpose of the off-diagonal block
  for(iint n=0;n<mesh->Np;++n){
    for(iint m=0;m<mesh->Np;++m){
      int id  = n*patchNp + mesh->Np + m;
      int idT = mesh->Np*patchNp + m*patchNp + n;

      A[idT] = A[id];
    }
  }
}


