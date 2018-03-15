#include "ellipticTri2D.h"

void matrixInverse(int N, dfloat *A);

dfloat matrixConditionNumber(int N, dfloat *A);

void BuildFaceIpdgPatchAx(solver_t *solver, mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong face, dfloat *A);

void BuildFaceBRdgPatchAx(solver_t *solver, mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong face, dfloat *A);


void ellipticBuildFacePatchesTri2D(solver_t *solver, mesh2D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   dlong *Npatches, dlong **patchesIndex, dfloat **patchesInvA,
                                   const char *options){

  if(!basis) { // default to degree N Lagrange basis
    basisNp = mesh->Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(int n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

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

  //We need the halo element's EToB flags to make the patch matrices
  if (mesh->totalHaloPairs) {
    mesh->EToB = (int *) realloc(mesh->EToB,(mesh->Nelements+mesh->totalHaloPairs)*sizeof(int));
    int *idSendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, sizeof(int), mesh->EToB, idSendBuffer, mesh->EToB + mesh->Nelements);
    free(idSendBuffer);
  }

  int NpatchElements = 2;
  int patchNp = mesh->Np*NpatchElements;

  //build a list of all face pairs
  mesh->NfacePairs=0;
  for (dlong eM=0; eM<mesh->Nelements;eM++) {
    for (int f=0;f<mesh->Nfaces;f++) {
      dlong eP = mesh->EToE[eM*mesh->Nfaces+f];

      if (eM<eP) mesh->NfacePairs++;
    }
  }

  mesh->EToFPairs = (dlong *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfaces,sizeof(dlong));
  mesh->FPairsToE = (dlong *) calloc(2*mesh->NfacePairs,sizeof(dlong));
  mesh->FPairsToF = (int *) calloc(2*mesh->NfacePairs,sizeof(int));

  //fill with -1
  for (dlong n=0;n<(mesh->Nelements+mesh->totalHaloPairs)*mesh->Nfaces;n++) {
    mesh->EToFPairs[n] = -1;
  }

  dlong cnt=0;
  for (dlong eM=0; eM<mesh->Nelements;eM++) {
    for (int fM=0;fM<mesh->Nfaces;fM++) {
      dlong eP = mesh->EToE[eM*mesh->Nfaces+fM];

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

  refMesh->EToV = (hlong*) calloc(NpatchElements*mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0*mesh->Nverts+0] = 0;
  refMesh->EToV[0*mesh->Nverts+1] = 1;
  refMesh->EToV[0*mesh->Nverts+2] = 2;

  refMesh->EToV[1*mesh->Nverts+0] = 1;
  refMesh->EToV[1*mesh->Nverts+1] = 0;
  refMesh->EToV[1*mesh->Nverts+2] = 3;

  refMesh->EToB = (int*) calloc(NpatchElements*mesh->Nfaces,sizeof(int));
  for (int n=0;n<NpatchElements*mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTri2D(refMesh, mesh->N);
  meshPhysicalNodesTri2D(refMesh);
  meshGeometricFactorsTri2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsTri2D(refMesh);

  //build a list of all face pairs
  refMesh->NfacePairs=1;

  refMesh->EToFPairs = (dlong *) calloc(2*mesh->Nfaces,sizeof(dlong));
  refMesh->FPairsToE = (dlong *) calloc(2,sizeof(dlong));
  refMesh->FPairsToF = (int *)  calloc(2,sizeof(int));

  //fill with -1
  for (int n=0;n<2*mesh->Nfaces;n++)  refMesh->EToFPairs[n] = -1;

  refMesh->FPairsToE[0] = 0;
  refMesh->FPairsToE[1] = 1;
  refMesh->FPairsToF[0] = 0;
  refMesh->FPairsToF[1] = 0;
  refMesh->EToFPairs[0] = 0;
  refMesh->EToFPairs[refMesh->Nfaces] = 1;



  int Nperm = pow(mesh->Nfaces,2);//all possible configuration of neighbours
  (*Npatches) = Nperm;
  dlong refPatches = 0;

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(Nperm*patchNp*patchNp, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->NfacePairs, sizeof(dlong));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));


  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  if (strstr(options,"IPDG")) {
    BuildFaceIpdgPatchAx(solver, refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  } else if (strstr(options,"BRDG")) {
    BuildFaceBRdgPatchAx(solver, refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  }
#if 0
  for (int n=0;n<mesh->Np*2;n++) {
    for (int m=0;m<mesh->Np*2;m++) {
      printf("%4.2f \t", refPatchInvA[m+n*mesh->Np*2]);
    }
    printf("\n");
  }
  printf("\n");
#endif
  matrixInverse(patchNp, refPatchInvA);

  //store the permutations of the reference patch
  int blkCounter = 1;
  int *permIndex = (int*) calloc(patchNp, sizeof(int));
  for(int blk=1;blk<Nperm;++blk){
    int f0 = blk%mesh->Nfaces;
    int f1 = (blk/mesh->Nfaces)%mesh->Nfaces;

    int r0 = (f0==0) ? 0: mesh->Nfaces-f0;
    int r1 = (f1==0) ? 0: mesh->Nfaces-f1;

    for(int n=0;n<mesh->Np;++n){
      permIndex[n+0*mesh->Np] = 0*mesh->Np + mesh->rmapP[r0*mesh->Np+n];
      permIndex[n+1*mesh->Np] = 1*mesh->Np + mesh->rmapP[r1*mesh->Np+n];
    }

    for(int n=0;n<patchNp;++n){
      for(int m=0;m<patchNp;++m){
        int pn = permIndex[n];
        int pm = permIndex[m];
        (*patchesInvA)[blkCounter*patchNp*patchNp + n*patchNp + m] = refPatchInvA[pn*patchNp+pm]; // maybe need to switch map
      }
    }
    ++blkCounter;
  }

  // loop over all elements
  for(dlong face=0;face<mesh->NfacePairs;++face){

    //build the patch A matrix for this element
    if (strstr(options,"IPDG")) {
      BuildFaceIpdgPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, face, patchA);
    } else if (strstr(options,"BRDG")) {
      BuildFaceBRdgPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, face, patchA);
    }
#if 0
    for (int n=0;n<mesh->Np*2;n++) {
      for (int m=0;m<mesh->Np*2;m++) {
        printf("%4.2f \t", patchA[m+n*mesh->Np*2]);
      }
      printf("\n");
    }
    printf("\n");
#endif
    dlong eM = mesh->FPairsToE[2*face+0];
    dlong eP = mesh->FPairsToE[2*face+1];
    int fM = mesh->FPairsToF[2*face+0];
    int fP = mesh->FPairsToF[2*face+1];

    if (eP >=0) {
      dlong eM0 = mesh->EToE[eM*mesh->Nfaces+0];
      dlong eM1 = mesh->EToE[eM*mesh->Nfaces+1];
      dlong eM2 = mesh->EToE[eM*mesh->Nfaces+2];

      dlong eP0 = mesh->EToE[eP*mesh->Nfaces+0];
      dlong eP1 = mesh->EToE[eP*mesh->Nfaces+1];
      dlong eP2 = mesh->EToE[eP*mesh->Nfaces+2];

      if(eM0>=0 && eM1>=0 && eM2>=0 &&
          eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interiour patch
        int blk = fM + mesh->Nfaces*fP;

        refPatchInvA = *patchesInvA + blk*patchNp*patchNp;

        //hit the patch with the reference inverse
        for(int n=0;n<patchNp;++n){
          for(int m=0;m<patchNp;++m){
            invRefAA[n*patchNp+m] = 0.;
            for (int k=0;k<patchNp;k++) {
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
    for(int n=0;n<patchNp;++n){
      for(int m=0;m<patchNp;++m){
        dlong id = ((*Npatches)-1)*patchNp*patchNp + n*patchNp + m;
        (*patchesInvA)[id] = patchA[n*patchNp+m];
      }
    }

    (*patchesIndex)[face] = (*Npatches)-1;
  }

  printf("saving "dlongFormat" full patches\n",*Npatches);
  printf("using "dlongFormat" reference patches\n", refPatches);

  free(patchA); free(invRefAA);
  free(MS);
}

void BuildFaceIpdgPatchAx(solver_t *solver, mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong face, dfloat *A) {

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
  for(int N=0;N<NpatchElements;++N){
    //element number
    dlong e = mesh->FPairsToE[2*face+N];

    dlong vbase = e*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
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

    for (int fM=0;fM<mesh->Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(e*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];

      int bc = mesh->EToB[fM+mesh->Nfaces*e]; //raw boundary flag

      dfloat penalty = tau*hinv;

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

      // mass matrix for this face
      dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<mesh->Nfp;++n){
        for(int m=0;m<mesh->Nfp;++m){
          int nM = mesh->faceNodes[fM*mesh->Nfp+n];
          int mM = mesh->faceNodes[fM*mesh->Nfp+m];
          int id = N*mesh->Np*patchNp + nM*patchNp + N*mesh->Np + mM;

          // OP11 = OP11 + 0.5*( gtau*mmE )
          dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
          A[id] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<mesh->Nfp;++n){
        for(int m=0;m<mesh->Np;++m){
          int nM = mesh->faceNodes[fM*mesh->Nfp+n];

          for(int i=0;i<mesh->Nfp;++i){
            int iM = mesh->faceNodes[fM*mesh->Nfp+i];

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

      for(int n=0;n<mesh->Np;++n){
        for(int m=0;m<mesh->Nfp;++m){
          int mM = mesh->faceNodes[fM*mesh->Nfp+m];

          for(int i=0;i<mesh->Nfp;++i){
            int iM = mesh->faceNodes[fM*mesh->Nfp+i];

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

    //add the rank boost for the allNeumann Poisson problem
    if (solver->allNeumann) {
      for(int n=0;n<mesh->Np;++n){
        for(int m=0;m<mesh->Np;++m){ 
          int id = N*mesh->Np*patchNp + n*patchNp + N*mesh->Np + m;
          A[id] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
        }
      }
    }
  }

  //now the off-diagonal
  dlong eM = mesh->FPairsToE[2*face+0];
  dlong eP = mesh->FPairsToE[2*face+1];
  int fM = mesh->FPairsToF[2*face+0];

  dlong sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
  dfloat nx = mesh->sgeo[sid+NXID];
  dfloat ny = mesh->sgeo[sid+NYID];
  dfloat sJ = mesh->sgeo[sid+SJID];
  dfloat hinv = mesh->sgeo[sid+IHID];

  dlong vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];

  dlong vbaseP = eP*mesh->Nvgeo;
  dfloat drdxP = mesh->vgeo[vbaseP+RXID];
  dfloat drdyP = mesh->vgeo[vbaseP+RYID];
  dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
  dfloat dsdyP = mesh->vgeo[vbaseP+SYID];

  dfloat penalty = tau*hinv;

  // mass matrix for this face
  dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

  // penalty term just involves face nodes
  for(int n=0;n<mesh->Nfp;++n){
    for(int m=0;m<mesh->Nfp;++m){
      int nM = mesh->faceNodes[fM*mesh->Nfp+n];
      // int mM = mesh->faceNodes[fM*mesh->Nfp+m];

      dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];

      // neighbor penalty term
      dlong idM = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+m;
      int mP = (int) (mesh->vmapP[idM]%mesh->Np);

      int id = nM*patchNp + mesh->Np + mP;

      // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
      A[id] += -0.5*penalty*MSfnm;
    }
  }

  // now add differential surface terms
  for(int n=0;n<mesh->Nfp;++n){
    for(int m=0;m<mesh->Np;++m){
      int nM = mesh->faceNodes[fM*mesh->Nfp+n];

      for(int i=0;i<mesh->Nfp;++i){
        // int iM = mesh->faceNodes[fM*mesh->Nfp+i];
        int iP = (int) (mesh->vmapP[i + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np);

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

  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Nfp;++m){
      // int mM = mesh->faceNodes[fM*mesh->Nfp+m];
      int mP = (int) (mesh->vmapP[m + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np);

      for(int i=0;i<mesh->Nfp;++i){
        int iM = mesh->faceNodes[fM*mesh->Nfp+i];

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

  if (solver->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ 
        int id = n*patchNp + mesh->Np + m;
        A[id] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
      }
    }
  }

  //write the transpose of the off-diagonal block
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      int id  = n*patchNp + mesh->Np + m;
      int idT = mesh->Np*patchNp + m*patchNp + n;

      A[idT] = A[id];
    }
  }
}



void BuildFaceBRdgPatchAx(solver_t *solver, mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong face, dfloat *A) {

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;

  int NpatchElements = 2;
  int patchNp = NpatchElements*Np;

  // Extract patches
  // B  a
  // a' C

  //zero out the matrix
  for (int n=0;n<patchNp*patchNp;n++) A[n] = 0.;

  int GblockSize = Np*Np*(Nfaces+1);

  /* Construct gradient as block matrix */
  dfloat  *Gx = (dfloat *) calloc(GblockSize*NpatchElements,sizeof(dfloat));
  dfloat  *Gy = (dfloat *) calloc(GblockSize*NpatchElements,sizeof(dfloat));

  //start with diagonals
  for(int N=0;N<NpatchElements;++N){
    //element number
    dlong e = mesh->FPairsToE[2*face+N];

    dlong vbase = e*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];

    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        Gx[m+n*Np+N*GblockSize] = drdx*mesh->Dr[m+n*Np]+dsdx*mesh->Ds[m+n*Np];
        Gy[m+n*Np+N*GblockSize] = drdy*mesh->Dr[m+n*Np]+dsdy*mesh->Ds[m+n*Np];
      }
    }

    for (int fM=0;fM<Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(e*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat invJ = mesh->sgeo[sid+IJID];

      dlong eP = mesh->EToE[e*Nfaces+fM];
      int fP = mesh->EToF[e*Nfaces+fM];
      dfloat sw = 1.f; //guard against unconnected elements (happens in reference patch)
      if (eP < 0) {eP = e; sw = 0;}
      if (fP < 0) fP = fM;

      // load surface geofactors for neighbor's face
      dlong sidP = mesh->Nsgeo*(eP*Nfaces+fP);
      dfloat nxP = mesh->sgeo[sidP+NXID];
      dfloat nyP = mesh->sgeo[sidP+NYID];
      dfloat sJP = mesh->sgeo[sidP+SJID];
      dfloat invJP = mesh->sgeo[sidP+IJID];

      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+Nfaces*e]; //raw boundary flag
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

      // lift term
      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh->faceNodes[fM*Nfp+m];

          dlong idP = eP*Nfp*Nfaces+fP*Nfp+m;          
          int mP = (int) (mesh->vmapP[idP]%Np);

          dfloat LIFTfnmM = sJ*invJ*mesh->LIFT[m + fM*Nfp + n*Nfp*Nfaces];
          dfloat LIFTfnmP = sJP*invJP*mesh->LIFT[m + fP*Nfp + n*Nfp*Nfaces];

          // G = sJ/J*LIFT*n*[[ uP-uM ]]
          Gx[mM+n*Np+N*GblockSize] += -0.5*(1-bcN)*(1+bcD)*nx*LIFTfnmM;
          Gy[mM+n*Np+N*GblockSize] += -0.5*(1-bcN)*(1+bcD)*ny*LIFTfnmM;

          Gx[mP+n*Np+(fM+1)*Np*Np+N*GblockSize] +=  0.5*sw*(1-bcN)*(1-bcD)*nxP*LIFTfnmP;
          Gy[mP+n*Np+(fM+1)*Np*Np+N*GblockSize] +=  0.5*sw*(1-bcN)*(1-bcD)*nyP*LIFTfnmP;
        }
      }
    }
  }

  //start with diagonals
  for(int N=0;N<NpatchElements;++N){
    //element number
    dlong e = mesh->FPairsToE[2*face+N];

    dlong vbase = e*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        int id = N*Np*patchNp + n*patchNp + N*Np + m;

        A[id]  = J*lambda*mesh->MM[n*Np+m];
        A[id] += J*drdx*drdx*mesh->Srr[n*Np+m];
        A[id] += J*drdx*dsdx*mesh->Srs[n*Np+m];
        A[id] += J*dsdx*drdx*mesh->Ssr[n*Np+m];
        A[id] += J*dsdx*dsdx*mesh->Sss[n*Np+m];

        A[id] += J*drdy*drdy*mesh->Srr[n*Np+m];
        A[id] += J*drdy*dsdy*mesh->Srs[n*Np+m];
        A[id] += J*dsdy*drdy*mesh->Ssr[n*Np+m];
        A[id] += J*dsdy*dsdy*mesh->Sss[n*Np+m];
      }
    }

    for (int fM=0;fM<Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(e*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];

      int bc = mesh->EToB[fM+Nfaces*e]; //raw boundary flag
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

      // mass matrix for this face
      dfloat *MSf = MS+fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Nfp;++m){
          int nM = mesh->faceNodes[fM*Nfp+n];
          int mM = mesh->faceNodes[fM*Nfp+m];
          int id = N*Np*patchNp + nM*patchNp + N*Np + mM;

          // OP11 = OP11 + 0.5*( gtau*mmE )
          dfloat MSfnm = sJ*MSf[n*Nfp+m];
          A[id] += 0.5*(1.-bcN)*(1.+bcD)*tau*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Np;++m){
          int nM = mesh->faceNodes[fM*Nfp+n];

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];
            int iP = (int) (mesh->vmapP[i+fM*Nfp+e*Nfp*Nfaces]%Np);

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            dfloat DxMim = Gx[m+iM*Np+N*GblockSize];
            dfloat DyMim = Gy[m+iM*Np+N*GblockSize];

            dfloat DxPim = Gx[m+iP*Np+(fM+1)*Np*Np+N*GblockSize];
            dfloat DyPim = Gy[m+iP*Np+(fM+1)*Np*Np+N*GblockSize];

            int id = N*Np*patchNp + nM*patchNp + N*Np + m;

            A[id] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            A[id] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

            A[id] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
            A[id] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
          }
        }
      }

      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh->faceNodes[fM*Nfp+m];

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

            int id = N*Np*patchNp + n*patchNp + N*Np + mM;

            // OP11 = OP11 + (- Dn1'*mmE );
            A[id] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            A[id] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (solver->allNeumann) {
      for(int n=0;n<Np;++n){
        for(int m=0;m<Np;++m){ 
          int id = N*Np*patchNp + n*patchNp + N*Np + m;
          A[id] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
        }
      }
    }
  }

  //now the off-diagonal
  dlong eM = mesh->FPairsToE[2*face+0];
  // dlong eP = mesh->FPairsToE[2*face+1];
  int fM = mesh->FPairsToF[2*face+0];
  int fP = mesh->FPairsToF[2*face+1];

  dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
  dfloat nx = mesh->sgeo[sid+NXID];
  dfloat ny = mesh->sgeo[sid+NYID];
  dfloat sJ = mesh->sgeo[sid+SJID];

  dlong vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];

  // mass matrix for this face
  dfloat *MSf = MS+fM*Nfp*Nfp;

  // penalty term just involves face nodes
  for(int n=0;n<Nfp;++n){
    for(int m=0;m<Nfp;++m){
      int nM = mesh->faceNodes[fM*Nfp+n];
      // int mM = mesh->faceNodes[fM*Nfp+m];

      dfloat MSfnm = sJ*MSf[n*Nfp+m];

      // neighbor penalty term
      dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
      int mP = (int) (mesh->vmapP[idM]%Np);

      int id = nM*patchNp + Np + mP;

      // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
      A[id] += -0.5*tau*MSfnm;
    }
  }

  // now add differential surface terms
  for(int n=0;n<Nfp;++n){
    for(int m=0;m<Np;++m){
      int nM = mesh->faceNodes[fM*Nfp+n];

      for(int i=0;i<Nfp;++i){
        int iM = mesh->faceNodes[fM*Nfp+i];
        int iP = (int) (mesh->vmapP[i + fM*Nfp+eM*Nfp*Nfaces]%Np);

        dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

        dfloat DxMim = Gx[m+iM*Np+(fP+1)*Np*Np+1*GblockSize];
        dfloat DyMim = Gy[m+iM*Np+(fP+1)*Np*Np+1*GblockSize];

        dfloat DxPim = Gx[m+iP*Np+1*GblockSize];
        dfloat DyPim = Gy[m+iP*Np+1*GblockSize];

        int id = nM*patchNp + Np + m;

        //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
        A[id] += -0.5*nx*MSfni*DxMim;
        A[id] += -0.5*ny*MSfni*DyMim;

        A[id] += -0.5*nx*MSfni*DxPim;
        A[id] += -0.5*ny*MSfni*DyPim;
      }
    }
  }

  for(int n=0;n<Np;++n){
    for(int m=0;m<Nfp;++m){
      // int mM = mesh->faceNodes[fM*Nfp+m];
      int mP = (int) (mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

      for(int i=0;i<Nfp;++i){
        int iM = mesh->faceNodes[fM*Nfp+i];

        dfloat MSfim = sJ*MSf[i*Nfp+m];

        dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
        dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

        int id = n*patchNp + Np + mP;

        //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
        A[id] +=  +0.5*nx*DxMin*MSfim;
        A[id] +=  +0.5*ny*DyMin*MSfim;
      }
    }
  }

  if (solver->allNeumann) {
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){ 
        int id = n*patchNp + Np + m;
        A[id] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
      }
    }
  }

  //write the transpose of the off-diagonal block
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      int id  = n*patchNp + Np + m;
      int idT = Np*patchNp + m*patchNp + n;

      A[idT] = A[id];
    }
  }

  free(Gx); free(Gy);
}


