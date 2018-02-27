#include "ellipticTri2D.h"

void matrixInverse(int N, dfloat *A);

dfloat matrixConditionNumber(int N, dfloat *A);

void BuildFullPatchAx(solver_t *solver, mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A);

void ellipticBuildFullPatchesTri2D(solver_t *solver, mesh2D* mesh, int basisNp, dfloat *basis,
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

  int NpatchElements = mesh->Nfaces+1;
  int patchNp = mesh->Np*NpatchElements;


  //build a mini mesh struct for the reference patch
  mesh2D *refMesh = (mesh2D*) calloc(1,sizeof(mesh2D));
  memcpy(refMesh,mesh,sizeof(mesh2D));

  //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x =        0., V4x =        0., V5x = 2., V6x = -2.;
  dfloat V1y = 0., V2y = 0., V3y =  sqrt(3.), V4y = -sqrt(3.), V5y = sqrt(3.), V6y =  sqrt(3.);

  refMesh->Nelements = NpatchElements;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));

  refMesh->EX[0*mesh->Nverts+0] = V1x;  refMesh->EY[0*mesh->Nverts+0] = V1y;
  refMesh->EX[0*mesh->Nverts+1] = V2x;  refMesh->EY[0*mesh->Nverts+1] = V2y;
  refMesh->EX[0*mesh->Nverts+2] = V3x;  refMesh->EY[0*mesh->Nverts+2] = V3y;

  refMesh->EX[1*mesh->Nverts+0] = V2x;  refMesh->EY[1*mesh->Nverts+0] = V2y;
  refMesh->EX[1*mesh->Nverts+1] = V1x;  refMesh->EY[1*mesh->Nverts+1] = V1y;
  refMesh->EX[1*mesh->Nverts+2] = V4x;  refMesh->EY[1*mesh->Nverts+2] = V4y;

  refMesh->EX[2*mesh->Nverts+0] = V3x;  refMesh->EY[2*mesh->Nverts+0] = V3y;
  refMesh->EX[2*mesh->Nverts+1] = V2x;  refMesh->EY[2*mesh->Nverts+1] = V2y;
  refMesh->EX[2*mesh->Nverts+2] = V5x;  refMesh->EY[2*mesh->Nverts+2] = V5y;

  refMesh->EX[3*mesh->Nverts+0] = V1x;  refMesh->EY[3*mesh->Nverts+0] = V1y;
  refMesh->EX[3*mesh->Nverts+1] = V3x;  refMesh->EY[3*mesh->Nverts+1] = V3y;
  refMesh->EX[3*mesh->Nverts+2] = V6x;  refMesh->EY[3*mesh->Nverts+2] = V6y;

  refMesh->EToV = (hlong*) calloc(NpatchElements*mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0*mesh->Nverts+0] = 0;
  refMesh->EToV[0*mesh->Nverts+1] = 1;
  refMesh->EToV[0*mesh->Nverts+2] = 2;

  refMesh->EToV[1*mesh->Nverts+0] = 1;
  refMesh->EToV[1*mesh->Nverts+1] = 0;
  refMesh->EToV[1*mesh->Nverts+2] = 3;

  refMesh->EToV[2*mesh->Nverts+0] = 2;
  refMesh->EToV[2*mesh->Nverts+1] = 1;
  refMesh->EToV[2*mesh->Nverts+2] = 4;

  refMesh->EToV[3*mesh->Nverts+0] = 0;
  refMesh->EToV[3*mesh->Nverts+1] = 2;
  refMesh->EToV[3*mesh->Nverts+2] = 5;

  refMesh->EToB = (int*) calloc(NpatchElements*mesh->Nfaces,sizeof(int));
  for (int n=0;n<NpatchElements*mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTri2D(refMesh, mesh->N);
  meshPhysicalNodesTri2D(refMesh);
  meshGeometricFactorsTri2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsTri2D(refMesh);

  int Nperm = pow(mesh->Nfaces,mesh->Nfaces);//all possible configureation of neighbours
  (*Npatches) = (dlong) Nperm;
  dlong refPatches = 0;

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(Nperm*patchNp*patchNp, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));


  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildFullPatchAx(solver, refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  matrixInverse(patchNp, refPatchInvA);

  //store the permutations of the reference patch
  int blkCounter = 1;
  int *permIndex = (int*) calloc(patchNp, sizeof(int));
  for(int blk=1;blk<Nperm;++blk){
    int f0 = blk%mesh->Nfaces;
    int f1 = (blk/mesh->Nfaces)%mesh->Nfaces;
    int f2= (blk/(mesh->Nfaces*mesh->Nfaces));

    int r0 = (f0==0) ? 0: mesh->Nfaces-f0;
    int r1 = (f1==0) ? 0: mesh->Nfaces-f1;
    int r2 = (f2==0) ? 0: mesh->Nfaces-f2;

    for(int n=0;n<mesh->Np;++n){
      permIndex[n+0*mesh->Np] = 0*mesh->Np + n;
      permIndex[n+1*mesh->Np] = 1*mesh->Np + mesh->rmapP[r0*mesh->Np+n];
      permIndex[n+2*mesh->Np] = 2*mesh->Np + mesh->rmapP[r1*mesh->Np+n];
      permIndex[n+3*mesh->Np] = 3*mesh->Np + mesh->rmapP[r2*mesh->Np+n];
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
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildFullPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, eM, patchA);

    dlong eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    dlong eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    dlong eP2 = mesh->EToE[eM*mesh->Nfaces+2];

    int fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    int fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    int fP2 = mesh->EToF[eM*mesh->Nfaces+2];

    if(eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interiour patch
      int blk = fP0 + mesh->Nfaces*fP1 + mesh->Nfaces*mesh->Nfaces*fP2;

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

      // printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < rateTolerance) {
        (*patchesIndex)[eM] = blk;
        refPatches++;
        continue;
      }
    }

    ++(*Npatches);
    *patchesInvA = (dfloat*) realloc(*patchesInvA, (*Npatches)*patchNp*patchNp*sizeof(dfloat));

    matrixInverse(patchNp, patchA);

    //copy inverse into patchesInvA
    for(int n=0;n<patchNp;++n){
      for(int m=0;m<patchNp;++m){
        int id = ((*Npatches)-1)*patchNp*patchNp + n*patchNp + m;
        (*patchesInvA)[id] = patchA[n*patchNp+m];
      }
    }

    (*patchesIndex)[eM] = (*Npatches)-1;
  }

  printf("saving %d full patches\n",*Npatches);
  printf("using %d reference patches\n", refPatches);

  free(refMesh);

  free(patchA); free(invRefAA);
  free(MS);
}


void BuildFullPatchAx(solver_t *solver, mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A) {

  int NpatchElements = mesh->Nfaces+1;
  int patchNp = NpatchElements*mesh->Np;

  // Extract patches
    // *  a b c
    // a' * 0 0
    // b' 0 * 0
    // c' 0 0 *

  //zero out the matrix
  for (int n=0;n<patchNp*patchNp;n++) A[n] = 0.;

  //make sure the diagonal is at least identity
  for (int n=0;n<patchNp;n++) A[n*patchNp+n] = 1.;

  //start with diagonals
  for(int N=0;N<NpatchElements;++N){
    dlong e = (N==0) ? eM : mesh->EToE[mesh->Nfaces*eM+N-1];

    if (e<0) continue; //skip this block if this is a boundary face

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

    //add the rank boost for the allNeumann Poisson problem
    if (solver->allNeumann) {
      for(int n=0;n<mesh->Np;++n){
        for(int m=0;m<mesh->Np;++m){ 
          A[n*mesh->Np+m] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
        }
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
  }

  //now the off-diagonals
  for (int fM=0;fM<mesh->Nfaces;fM++) {

    dlong eP = mesh->EToE[eM*mesh->Nfaces+fM];
    if (eP < 0) continue; //skip this block if this is a boundary face

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
    dfloat J = mesh->vgeo[vbase+JID];

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
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];

        dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];

        // neighbor penalty term
        dlong idM = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+m;
        int mP = (int) (mesh->vmapP[idM]%mesh->Np);

        int id = nM*patchNp + (fM+1)*mesh->Np + mP;

        // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
        A[id] += -0.5*penalty*MSfnm;
      }
    }

    // now add differential surface terms
    for(int n=0;n<mesh->Nfp;++n){
      for(int m=0;m<mesh->Np;++m){
        int nM = mesh->faceNodes[fM*mesh->Nfp+n];

        for(int i=0;i<mesh->Nfp;++i){
          int iM = mesh->faceNodes[fM*mesh->Nfp+i];
          int iP = (int) (mesh->vmapP[i + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np);

          dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in

          dfloat DxPim = drdxP*mesh->Dr[iP*mesh->Np+m] + dsdxP*mesh->Ds[iP*mesh->Np+m];
          dfloat DyPim = drdyP*mesh->Dr[iP*mesh->Np+m] + dsdyP*mesh->Ds[iP*mesh->Np+m];

          int id = nM*patchNp + (fM+1)*mesh->Np + m;

          //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
          A[id] += -0.5*nx*MSfni*DxPim;
          A[id] += -0.5*ny*MSfni*DyPim;
        }
      }
    }

    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Nfp;++m){
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];
        int mP = (int) (mesh->vmapP[m + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np);

        for(int i=0;i<mesh->Nfp;++i){
          int iM = mesh->faceNodes[fM*mesh->Nfp+i];

          dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

          dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
          dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];

          int id = n*patchNp + (fM+1)*mesh->Np + mP;

          //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
          A[id] +=  +0.5*nx*DxMin*MSfim;
          A[id] +=  +0.5*ny*DyMin*MSfim;
        }
      }
    }

    //write the transpose of the off-diagonal block
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        int id  = n*patchNp + (fM+1)*mesh->Np + m;
        int idT = (fM+1)*mesh->Np*patchNp + m*patchNp + n;

        A[idT] = A[id];
      }
    }
  }
}


void matrixInverse(int N, dfloat *A){
  int lwork = N*N;
  int info;

  // compute inverse mass matrix
  double *tmpInvA = (double*) calloc(N*N, sizeof(double));

  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(int n=0;n<N*N;++n){
    tmpInvA[n] = A[n];
  }

  dgetrf_ (&N, &N, tmpInvA, &N, ipiv, &info);
  dgetri_ (&N, tmpInvA, &N, ipiv, work, &lwork, &info);

  if(info)
    printf("inv: dgetrf/dgetri reports info = %d when inverting matrix\n", info);

  for(int n=0;n<N*N;++n)
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

  int *ipiv = (int*) calloc(N, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));
  int  *iwork = (int*) calloc(N, sizeof(int));

  for(int n=0;n<N*N;++n){
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