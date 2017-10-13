#include "ellipticTet3D.h"

typedef struct {

  int face;
  int signature[4];
  iint id;

} refPatch_t;


void matrixInverse(int N, dfloat *A);

dfloat matrixConditionNumber(int N, dfloat *A);

void BuildFacePatchAx(solver_t *solver, mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint face, dfloat *A);

void BuildReferenceFacePatch(solver_t *solver, mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint face, int *signature, dfloat *A);

iint getFacePatchIndex(refPatch_t *referencePatchList, iint numRefPatches, int face, int *signature);




void ellipticBuildFacePatchesIpdgTet3D(solver_t *solver, mesh3D *mesh, iint basisNp, dfloat *basis,
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


  (*Npatches) = 0;
  int numRefPatches=0;
  int refPatches = 0;

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc((*Npatches)*patchNp*patchNp, sizeof(dfloat));
  *patchesIndex = (iint*) calloc(mesh->NfacePairs, sizeof(iint));

  refPatch_t *referencePatchList = (refPatch_t *) calloc(numRefPatches,sizeof(refPatch_t));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));

  dfloat *refPatchInvA;

  // loop over all face pairs
  for(iint face=0;face<mesh->NfacePairs;++face){

    //build the patch A matrix for this element
    BuildFacePatchAx(solver, mesh, basis, tau, lambda, BCType, MS, face, patchA);

    iint eM = mesh->FPairsToE[2*face+0];
    iint eP = mesh->FPairsToE[2*face+1];
    iint fM = mesh->FPairsToF[2*face+0];
    iint fP = mesh->FPairsToF[2*face+1];

    iint eM0 = mesh->EToE[eM*mesh->Nfaces+0];
    iint eM1 = mesh->EToE[eM*mesh->Nfaces+1];
    iint eM2 = mesh->EToE[eM*mesh->Nfaces+2];
    iint eM3 = mesh->EToE[eM*mesh->Nfaces+3];

    iint eP0 = mesh->EToE[eP*mesh->Nfaces+0];
    iint eP1 = mesh->EToE[eP*mesh->Nfaces+1];
    iint eP2 = mesh->EToE[eP*mesh->Nfaces+2];
    iint eP3 = mesh->EToE[eP*mesh->Nfaces+3];

    if(eM0>=0 && eM1>=0 && eM2>=0 && eM3>=0 &&
       eP0>=0 && eP1>=0 && eP2>=0 && eP3>=0){ //check if this is an interiour patch
      
      //get the vertices
      iint *vM = mesh->EToV + eM*mesh->Nverts;
      iint *vP = mesh->EToV + eP*mesh->Nverts;

      // intialize signature to -1
      int signature[4];
      for (int n=0;n<mesh->Nverts;n++) signature[n] = -1;

      for (int n=0;n<mesh->Nverts;n++) {
        for (int m=0;m<mesh->Nverts;m++) {
          if (vP[m] == vM[n]) signature[m] = n; 
        }
      }

      iint index = getFacePatchIndex(referencePatchList,numRefPatches,fM,signature);
      if (index<0) {      
        //build the reference patch for this signature
        ++(*Npatches);
        numRefPatches++;
        *patchesInvA = (dfloat*) realloc(*patchesInvA, (*Npatches)*patchNp*patchNp*sizeof(dfloat));
        referencePatchList = (refPatch_t *) realloc(referencePatchList,numRefPatches*sizeof(refPatch_t));

        referencePatchList[numRefPatches-1].face = fM; 
        referencePatchList[numRefPatches-1].id = (*Npatches)-1;
        for (int n=0;n<mesh->Nverts;n++) 
          referencePatchList[numRefPatches-1].signature[n] = signature[n];

        refPatchInvA = *patchesInvA + ((*Npatches)-1)*patchNp*patchNp;

        // printf("Building reference patch for the face %d, with signature [%d,%d,%d,%d] \n", fM, signature[0], signature[1], signature[2],signature[3]);

        BuildReferenceFacePatch(solver, mesh, basis, tau, lambda, BCType, MS, fM, signature, refPatchInvA); 
        matrixInverse(patchNp, refPatchInvA);        
        index = (*Npatches)-1;
      }

      refPatchInvA = *patchesInvA + index*patchNp*patchNp;

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

      printf("Face pair %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      //if the convergence rate is good, use the reference patch, and skip adding this patch
      if (rate < rateTolerance) {
        (*patchesIndex)[face] = index;
        refPatches++;
        continue;
      }
    }
    //add this patch to the patch list
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


void BuildFacePatchAx(solver_t *solver, mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
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
        int id = N*mesh->Np*patchNp + n*patchNp + N*mesh->Np + m;

        A[id]  = J*lambda*mesh->MM[n*mesh->Np+m];
        A[id] += J*G00*mesh->Srr[n*mesh->Np+m];
        A[id] += J*G01*mesh->Srs[n*mesh->Np+m];
        A[id] += J*G02*mesh->Srt[n*mesh->Np+m];
        A[id] += J*G10*mesh->Ssr[n*mesh->Np+m];
        A[id] += J*G11*mesh->Sss[n*mesh->Np+m];
        A[id] += J*G12*mesh->Sst[n*mesh->Np+m];
        A[id] += J*G20*mesh->Str[n*mesh->Np+m];
        A[id] += J*G21*mesh->Sts[n*mesh->Np+m];
        A[id] += J*G22*mesh->Stt[n*mesh->Np+m];
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (solver->allNeumann) {
      for(iint n=0;n<mesh->Np;++n){
        for(iint m=0;m<mesh->Np;++m){ 
          A[n*mesh->Np+m] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
        }
      }
    }

    for (iint fM=0;fM<mesh->Nfaces;fM++) {
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(e*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat nz = mesh->sgeo[sid+NZID];
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

            dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+m] + dsdx*mesh->Ds[iM*mesh->Np+m]+ dtdx*mesh->Dt[iM*mesh->Np+m];
            dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+m] + dsdy*mesh->Ds[iM*mesh->Np+m]+ dtdy*mesh->Dt[iM*mesh->Np+m];
            dfloat DzMim = drdz*mesh->Dr[iM*mesh->Np+m] + dsdz*mesh->Ds[iM*mesh->Np+m]+ dtdz*mesh->Dt[iM*mesh->Np+m];

            int id = N*mesh->Np*patchNp + nM*patchNp + N*mesh->Np + m;

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            A[id] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            A[id] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
            A[id] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;
          }
        }
      }

      for(iint n=0;n<mesh->Np;++n){
        for(iint m=0;m<mesh->Nfp;++m){
          iint mM = mesh->faceNodes[fM*mesh->Nfp+m];

          for(iint i=0;i<mesh->Nfp;++i){
            iint iM = mesh->faceNodes[fM*mesh->Nfp+i];

            dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n]+ dtdx*mesh->Dt[iM*mesh->Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n]+ dtdy*mesh->Dt[iM*mesh->Np+n];
            dfloat DzMin = drdz*mesh->Dr[iM*mesh->Np+n] + dsdz*mesh->Ds[iM*mesh->Np+n]+ dtdz*mesh->Dt[iM*mesh->Np+n];

            int id = N*mesh->Np*patchNp + n*patchNp + N*mesh->Np + mM;

            // OP11 = OP11 + (- Dn1'*mmE );
            A[id] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            A[id] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
            A[id] +=  -0.5*nz*(1+bcD)*(1-bcN)*DzMin*MSfim;
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
  dfloat nz = mesh->sgeo[sid+NZID];
  dfloat sJ = mesh->sgeo[sid+SJID];
  dfloat hinv = mesh->sgeo[sid+IHID];

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

  iint vbaseP = eP*mesh->Nvgeo;
  dfloat drdxP = mesh->vgeo[vbaseP+RXID];
  dfloat drdyP = mesh->vgeo[vbaseP+RYID];
  dfloat drdzP = mesh->vgeo[vbaseP+RZID];
  dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
  dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
  dfloat dsdzP = mesh->vgeo[vbaseP+SZID];
  dfloat dtdxP = mesh->vgeo[vbaseP+TXID];
  dfloat dtdyP = mesh->vgeo[vbaseP+TYID];
  dfloat dtdzP = mesh->vgeo[vbaseP+TZID];

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

        dfloat DxPim = drdxP*mesh->Dr[iP*mesh->Np+m] + dsdxP*mesh->Ds[iP*mesh->Np+m]+ dtdxP*mesh->Dt[iP*mesh->Np+m];
        dfloat DyPim = drdyP*mesh->Dr[iP*mesh->Np+m] + dsdyP*mesh->Ds[iP*mesh->Np+m]+ dtdyP*mesh->Dt[iP*mesh->Np+m];
        dfloat DzPim = drdzP*mesh->Dr[iP*mesh->Np+m] + dsdzP*mesh->Ds[iP*mesh->Np+m]+ dtdzP*mesh->Dt[iP*mesh->Np+m];

        int id = nM*patchNp + mesh->Np + m;

        //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
        A[id] += -0.5*nx*MSfni*DxPim;
        A[id] += -0.5*ny*MSfni*DyPim;
        A[id] += -0.5*nz*MSfni*DzPim;
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

        dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n]+ dtdx*mesh->Dt[iM*mesh->Np+n];
        dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n]+ dtdy*mesh->Dt[iM*mesh->Np+n];
        dfloat DzMin = drdz*mesh->Dr[iM*mesh->Np+n] + dsdz*mesh->Ds[iM*mesh->Np+n]+ dtdz*mesh->Dt[iM*mesh->Np+n];

        int id = n*patchNp + mesh->Np + mP;

        //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
        A[id] +=  +0.5*nx*DxMin*MSfim;
        A[id] +=  +0.5*ny*DyMin*MSfim;
        A[id] +=  +0.5*nz*DzMin*MSfim;
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

void BuildReferenceFacePatch(solver_t *solver, mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint face, int *signature, dfloat *A) {
  //build a mini mesh struct for the reference patch
  mesh3D *refMesh = (mesh3D*) calloc(1,sizeof(mesh3D));
  memcpy(refMesh,mesh,sizeof(mesh3D));

   //vertices of reference patch
  int Nv = 8;
  dfloat VX[8] = {-1, 1,      0,          0,           0,         5./3,        -5./3,         0};
  dfloat VY[8] = { 0, 0,sqrt(3.),  1/sqrt(3.),-7*sqrt(3.)/9, 8*sqrt(3.)/9, 8*sqrt(3.)/9, 1/sqrt(3.)};
  dfloat VZ[8] = { 0, 0,      0,2*sqrt(6.)/3, 4*sqrt(6.)/9, 4*sqrt(6.)/9, 4*sqrt(6.)/9,-2*sqrt(6.)/3};

  iint EToV[5*4] = {0,1,2,3,
                    0,2,1,7,
                    0,1,3,4,
                    1,2,3,5,
                    2,0,3,6};

  int NpatchElements = 2;                    
  refMesh->Nelements = NpatchElements;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));
  refMesh->EZ = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));

  refMesh->EToV = (iint*) calloc(NpatchElements*mesh->Nverts, sizeof(iint));


  for(int n=0;n<mesh->Nverts;++n){
    int v = EToV[n];
    refMesh->EX[n] = VX[v];
    refMesh->EY[n] = VY[v];
    refMesh->EZ[n] = VZ[v];
    refMesh->EToV[n] = v;
  } 

  for (int n=0;n<mesh->Nverts;n++) {
    if (signature[n]==-1) {
      //fill the missing vertex based ont he face number
      int v = EToV[(face+1)*mesh->Nverts+mesh->Nverts-1];  
      refMesh->EX[mesh->Nverts+n] = VX[v];
      refMesh->EY[mesh->Nverts+n] = VY[v];
      refMesh->EZ[mesh->Nverts+n] = VZ[v];
      refMesh->EToV[mesh->Nverts+n] = 4; //extra vert      
    } else {
      int v = signature[n];  
      refMesh->EX[mesh->Nverts+n] = VX[v];
      refMesh->EY[mesh->Nverts+n] = VY[v];
      refMesh->EZ[mesh->Nverts+n] = VZ[v];
      refMesh->EToV[mesh->Nverts+n] = v;      
    }
  }  

  refMesh->EToB = (iint*) calloc(NpatchElements*mesh->Nfaces,sizeof(iint));
  for (iint n=0;n<NpatchElements*mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  //build a list of all face pairs
  refMesh->NfacePairs=1;

  refMesh->EToFPairs = (iint *) calloc(2*mesh->Nfaces,sizeof(iint));
  refMesh->FPairsToE = (iint *) calloc(2,sizeof(iint));
  refMesh->FPairsToF = (int *)  calloc(2,sizeof(int));

  //fill with -1
  for (iint n=0;n<2*mesh->Nfaces;n++)  refMesh->EToFPairs[n] = -1;

  refMesh->FPairsToE[0] = 0;
  refMesh->FPairsToE[1] = 1;
  refMesh->FPairsToF[0] = face;
  refMesh->FPairsToF[1] = 0;
  refMesh->EToFPairs[0] = 0;
  refMesh->EToFPairs[refMesh->Nfaces] = 1;

  meshConnect(refMesh);
  meshLoadReferenceNodesTet3D(refMesh, mesh->N);
  meshPhysicalNodesTet3D(refMesh);
  meshGeometricFactorsTet3D(refMesh);
  meshConnectFaceNodes3D(refMesh);
  meshSurfaceGeometricFactorsTet3D(refMesh);

  //build this reference patch
  BuildFacePatchAx(solver, refMesh, basis, tau, lambda, BCType, MS, 0, A);

  free(refMesh->EX);
  free(refMesh->EY);
  free(refMesh->EZ);
  free(refMesh->EToV);

  free(refMesh);
}

iint getFacePatchIndex(refPatch_t *referencePatchList, iint numRefPatches, int face, int *signature) {

  iint index = -1;
  for (iint n=0;n<numRefPatches;n++) {
    if (referencePatchList[n].face == face) {
      if ((referencePatchList[n].signature[0] == signature[0]) &&
          (referencePatchList[n].signature[1] == signature[1]) && 
          (referencePatchList[n].signature[2] == signature[2]) &&
          (referencePatchList[n].signature[3] == signature[3])) {
        index = referencePatchList[n].id;
        break;
      }
    }
  }
  return index;
}





