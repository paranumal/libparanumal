#include "ellipticTet3D.h"

typedef struct {

  int signature[4*4];
  iint id;

} refPatch_t;

void matrixInverse(int N, dfloat *A);

dfloat matrixConditionNumber(int N, dfloat *A);


void BuildFullPatchAx(mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A);

void BuildReferenceFullPatch(mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, int *signature, dfloat *A);

iint getFullPatchIndex(refPatch_t *referencePatchList, iint numRefPatches, int *signature);



void ellipticBuildFullPatchesIpdgTet3D(mesh3D *mesh, iint basisNp, dfloat *basis,
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

  int NpatchElements = mesh->Nfaces+1;
  int patchNp = mesh->Np*NpatchElements;

  (*Npatches) = 0;
  int numRefPatches=0;
  int refPatches = 0;

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc((*Npatches)*patchNp*patchNp, sizeof(dfloat));
  *patchesIndex = (iint*) calloc(mesh->Nelements, sizeof(iint));

  refPatch_t *referencePatchList = (refPatch_t *) calloc(numRefPatches,sizeof(refPatch_t));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));

  dfloat *refPatchInvA;
  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildFullPatchAx(mesh, basis, tau, lambda, BCType, MS, eM, patchA);

    iint eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    iint eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    iint eP2 = mesh->EToE[eM*mesh->Nfaces+2];
    iint eP3 = mesh->EToE[eM*mesh->Nfaces+3];

    if(eP0>=0 && eP1>=0 && eP2>=0 && eP3>=0){ //check if this is an interiour patch

      //get the vertices
      iint *vM = mesh->EToV + eM*mesh->Nverts;
      iint *vP0 = mesh->EToV + eP0*mesh->Nverts;
      iint *vP1 = mesh->EToV + eP1*mesh->Nverts;
      iint *vP2 = mesh->EToV + eP2*mesh->Nverts;
      iint *vP3 = mesh->EToV + eP3*mesh->Nverts;

      iint *vP[4] = {vP0,vP1,vP2,vP3};

      // intialize signature to -1
      int signature[4*4];
      for (int n=0;n<mesh->Nfaces*mesh->Nverts;n++) signature[n] = -1;

      for (int f=0;f<mesh->Nfaces;f++) {
        for (int n=0;n<mesh->Nverts;n++) {
          for (int m=0;m<mesh->Nverts;m++) {
            if (vP[f][m] == vM[n]) signature[f*mesh->Nverts + m] = n; 
          }
        }
      }

      iint index = getFullPatchIndex(referencePatchList,numRefPatches,signature);
      if (index<0) {      
        //build the reference patch for this signature
        ++(*Npatches);
        numRefPatches++;
        *patchesInvA = (dfloat*) realloc(*patchesInvA, (*Npatches)*patchNp*patchNp*sizeof(dfloat));
        referencePatchList = (refPatch_t *) realloc(referencePatchList,numRefPatches*sizeof(refPatch_t)); 
        referencePatchList[numRefPatches-1].id = (*Npatches)-1;
        for (int n=0;n<mesh->Nverts*mesh->Nfaces;n++) 
          referencePatchList[numRefPatches-1].signature[n] = signature[n];

        refPatchInvA = *patchesInvA + ((*Npatches)-1)*patchNp*patchNp;

        // printf("Building reference patch with signature [%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d] \n", 
                  // signature[0], signature[1], signature[2],signature[3],
                  // signature[4], signature[5], signature[6],signature[7],
                  // signature[8], signature[9], signature[10],signature[11],
                  // signature[12], signature[13], signature[14],signature[15]);

        BuildReferenceFullPatch(mesh, basis, tau, lambda, BCType, MS, signature, refPatchInvA); 
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

      printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = index;
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

    (*patchesIndex)[eM] = (*Npatches)-1;
  }

  printf("saving %d full patches\n",*Npatches);
  printf("using %d reference patches\n", refPatches);

  free(patchA); free(invRefAA);
  free(MS);
}

void BuildFullPatchAx(mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A) {;

  int NpatchElements = mesh->Nfaces+1;
  int patchNp = NpatchElements*mesh->Np;

  // Extract patches
    // *  a b c d
    // a' * 0 0 0
    // b' 0 * 0 0
    // c' 0 0 * 0
    // d' 0 0 0 *

  //zero out the matrix
  for (int n=0;n<patchNp*patchNp;n++) A[n] = 0.;

  //make sure the diagonal is at least identity
  for (int n=0;n<patchNp;n++) A[n*patchNp+n] = 1.;

  //start with diagonals
  for(iint N=0;N<NpatchElements;++N){
    iint e = (N==0) ? eM : mesh->EToE[mesh->Nfaces*eM+N-1];

    if (e<0) continue; //skip this block if this is a boundary face

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

  //now the off-diagonals
  for (iint fM=0;fM<mesh->Nfaces;fM++) {

    iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
    if (eP < 0) continue; //skip this block if this is a boundary face

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

        int id = nM*patchNp + (fM+1)*mesh->Np + mP;

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

          int id = nM*patchNp + (fM+1)*mesh->Np + m;

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

          int id = n*patchNp + (fM+1)*mesh->Np + mP;

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
        int id  = n*patchNp + (fM+1)*mesh->Np + m;
        int idT = (fM+1)*mesh->Np*patchNp + m*patchNp + n;

        A[idT] = A[id];
      }
    }
  }
}

void BuildReferenceFullPatch(mesh3D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, int *signature, dfloat *A) {
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

  int NpatchElements = 5;                    
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

  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nverts;n++) {
      if (signature[f*mesh->Nverts + n]==-1) {
        //fill the missing vertex based ont he face number
        int v = EToV[(f+1)*mesh->Nverts+mesh->Nverts-1];  
        refMesh->EX[(f+1)*mesh->Nverts+n] = VX[v];
        refMesh->EY[(f+1)*mesh->Nverts+n] = VY[v];
        refMesh->EZ[(f+1)*mesh->Nverts+n] = VZ[v];
        refMesh->EToV[(f+1)*mesh->Nverts+n] = v; //extra vert      
      } else {
        int v = signature[f*mesh->Nverts + n];  
        refMesh->EX[(f+1)*mesh->Nverts+n] = VX[v];
        refMesh->EY[(f+1)*mesh->Nverts+n] = VY[v];
        refMesh->EZ[(f+1)*mesh->Nverts+n] = VZ[v];
        refMesh->EToV[(f+1)*mesh->Nverts+n] = v;      
      }
    }  
  }

  refMesh->EToB = (iint*) calloc(NpatchElements*mesh->Nfaces,sizeof(iint));
  for (iint n=0;n<NpatchElements*mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTet3D(refMesh, mesh->N);
  meshPhysicalNodesTet3D(refMesh);
  meshGeometricFactorsTet3D(refMesh);
  meshConnectFaceNodes3D(refMesh);
  meshSurfaceGeometricFactorsTet3D(refMesh);

  //build this reference patch
  BuildFullPatchAx(refMesh, basis, tau, lambda, BCType, MS, 0, A);
;
  free(refMesh->EX);
  free(refMesh->EY);
  free(refMesh->EZ);
  free(refMesh->EToV);
  free(refMesh->EToB);

  free(refMesh);
}

iint getFullPatchIndex(refPatch_t *referencePatchList, iint numRefPatches, int *signature) {

  iint index = -1;
  for (iint n=0;n<numRefPatches;n++) {
    bool match = true;
    for (int m=0;m<4*4;m++) {
      match = match && (referencePatchList[n].signature[n] == signature[n]);
    }
    if (match) {
      index = referencePatchList[n].id;
      break;
    }
  }
  return index;
}

void matrixInverse(int N, dfloat *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  double *tmpInvA = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (iint*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(iint n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  dgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  dgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: dgetrf/dgetri reports info = %d when inverting matrix\n", info);

  for(iint n=0;n<N*N;++n)
    A[n] = tmpInvA[n];

  free(work);
  free(ipiv);
  free(tmpInvA);
}

dfloat matrixConditionNumber(int N, dfloat *A) {

  int lwork = 4*N;
  int info;

  char norm = '1';

  double Acond;
  double Anorm;

  double *tmpLU = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (iint*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));
  int  *iwork = (int*) calloc(N, sizeof(int));

  for(iint n=0;n<N*N;++n){
    tmpLU[n] = (double) A[n];
  }

  //get the matrix norm of A
  Anorm = dlange_(&norm, &N, &N, tmpLU, &N, work);

  //compute LU factorization
  dgetrf_ (&N, &N, tmpLU, &N, ipiv, &info);

  //compute inverse condition number
  dgecon_(&norm, &N, tmpLU, &N, &Anorm, &Acond, work, iwork, &info);

  if(info)
    printf("inv: dgetrf/dgecon reports info = %d when computing condition number\n", info);

  free(work);
  free(iwork);
  free(ipiv);
  free(tmpLU);

  return (dfloat) 1.0/Acond;
}
