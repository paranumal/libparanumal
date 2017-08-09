#include "ellipticTri2D.h"

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

void BuildFullPatchAx(mesh2D *mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, int* refFaceMap, iint eM, dfloat *A) {

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
  for(iint N=0;N<NpatchElements;++N){
    iint e;

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
      //element number
      e = (N==0) ? eM : mesh->EToE[mesh->Nfaces*eM+N-1];

      if (e<0) continue; //skip this block if this is a boundary face

      vbase = e*mesh->Nvgeo;
      drdx = mesh->vgeo[vbase+RXID];
      drdy = mesh->vgeo[vbase+RYID];
      dsdx = mesh->vgeo[vbase+SXID];
      dsdy = mesh->vgeo[vbase+SYID];
      J = mesh->vgeo[vbase+JID];
    }

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
        sid = mesh->Nsgeo*(e*mesh->Nfaces+fM);
        nx = mesh->sgeo[sid+NXID];
        ny = mesh->sgeo[sid+NYID];
        sJ = mesh->sgeo[sid+SJID];
        hinv = mesh->sgeo[sid+IHID];

        bc = mesh->EToB[fM+mesh->Nfaces*e]; //raw boundary flag
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

  //now the off-diagonals
  for (iint fM=0;fM<mesh->Nfaces;fM++) {
    iint sid;
    dfloat nx;
    dfloat ny;
    dfloat sJ;
    dfloat hinv;

    iint vbase;
    dfloat drdx;
    dfloat drdy;
    dfloat dsdx;
    dfloat dsdy;
    dfloat J;

    iint eP = 0;
    iint vbaseP;
    dfloat drdxP;
    dfloat drdyP;
    dfloat dsdxP;
    dfloat dsdyP;

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

      drdx = 1.0;
      drdy = 1./sqrt(3.);
      dsdx = 0.;
      dsdy = 2./sqrt(3.);
      J = sqrt(3.)/2.;

      drdxP = 1.0;
      drdyP = 1./sqrt(3.);
      dsdxP = 0.;
      dsdyP = 2./sqrt(3.);
    } else {
      eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) continue; //skip this block if this is a boundary face

      sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
      nx = mesh->sgeo[sid+NXID];
      ny = mesh->sgeo[sid+NYID];
      sJ = mesh->sgeo[sid+SJID];
      hinv = mesh->sgeo[sid+IHID];

      vbase = eM*mesh->Nvgeo;
      drdx = mesh->vgeo[vbase+RXID];
      drdy = mesh->vgeo[vbase+RYID];
      dsdx = mesh->vgeo[vbase+SXID];
      dsdy = mesh->vgeo[vbase+SYID];
      J = mesh->vgeo[vbase+JID];

      vbaseP = eP*mesh->Nvgeo;
      drdxP = mesh->vgeo[vbaseP+RXID];
      drdyP = mesh->vgeo[vbaseP+RYID];
      dsdxP = mesh->vgeo[vbaseP+SXID];
      dsdyP = mesh->vgeo[vbaseP+SYID];

      // reset eP
      eP = mesh->EToE[eM*mesh->Nfaces+fM];
    }

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
        if(eP>=0){
          int mP;
          if(eM==-1) {
            mP = refFaceMap[fM*mesh->Nfp+m];
          } else {
            iint idM = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+m;
            mP = mesh->vmapP[idM]%mesh->Np;
          }

          int id = nM*patchNp + (fM+1)*mesh->Np + mP;

          // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
          A[id] += -0.5*penalty*MSfnm;
        }
      }
    }

    // now add differential surface terms
    for(iint n=0;n<mesh->Nfp;++n){
      for(iint m=0;m<mesh->Np;++m){
        int nM = mesh->faceNodes[fM*mesh->Nfp+n];

        for(iint i=0;i<mesh->Nfp;++i){
          int iM = mesh->faceNodes[fM*mesh->Nfp+i];
          int iP;

          if(eM==-1) {
            iP = refFaceMap[fM*mesh->Nfp+i];
          } else {
            iP = mesh->vmapP[i + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;
          }

          dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in

          if(eP>=0){
            dfloat DxPim = drdxP*mesh->Dr[iP*mesh->Np+m] + dsdxP*mesh->Ds[iP*mesh->Np+m];
            dfloat DyPim = drdyP*mesh->Dr[iP*mesh->Np+m] + dsdyP*mesh->Ds[iP*mesh->Np+m];

            int id = nM*patchNp + (fM+1)*mesh->Np + m;

            //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
            A[id] += -0.5*nx*MSfni*DxPim;
            A[id] += -0.5*ny*MSfni*DyPim;
          }
        }
      }
    }

    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Nfp;++m){
        int mM = mesh->faceNodes[fM*mesh->Nfp+m];
        int mP;

        if(eM==-1) {
          mP = refFaceMap[fM*mesh->Nfp+m];
        } else {
          mP = mesh->vmapP[m + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;
        }

        for(iint i=0;i<mesh->Nfp;++i){
          iint iM = mesh->faceNodes[fM*mesh->Nfp+i];

          dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];

          dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
          dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];

          if(eP>=0){
            int id = n*patchNp + (fM+1)*mesh->Np + mP;

            //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
            A[id] +=  +0.5*nx*DxMin*MSfim;
            A[id] +=  +0.5*ny*DyMin*MSfim;
          }
        }
      }
    }

    //write the transpose of the off-diagonal block
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Nfp;++m){
        int id  = n*patchNp + (fM+1)*mesh->Np + m;
        int idT = (fM+1)*mesh->Np*patchNp + m*patchNp + n;

        A[idT] = A[id];
      }
    }
  }
}

void ellipticBuildExactPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
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

  //We need the halo element's EToB flags to make the patch matrices
  if (mesh->totalHaloPairs) {
    mesh->EToB = (int *) realloc(mesh->EToB,(mesh->Nelements+mesh->totalHaloPairs)*sizeof(int));
    iint *idSendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, sizeof(int), mesh->EToB, idSendBuffer, mesh->EToB + mesh->Nelements);
    free(idSendBuffer);
  }

  int NpatchElements = mesh->Nfaces+1;
  int patchNp = mesh->Np*NpatchElements;

  //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x = -1., V4x =  1., V5x = 1., V6x = -3.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y = -3., V5y = 1., V6y =  1.;

  dfloat *EX = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));
  dfloat *EY = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));

  EX[0*mesh->Nverts+0] = V1x;  EY[0*mesh->Nverts+0] = V1y;
  EX[0*mesh->Nverts+1] = V2x;  EY[0*mesh->Nverts+1] = V2y;
  EX[0*mesh->Nverts+2] = V3x;  EY[0*mesh->Nverts+2] = V3y;

  EX[1*mesh->Nverts+0] = V2x;  EY[1*mesh->Nverts+0] = V2y;
  EX[1*mesh->Nverts+1] = V1x;  EY[1*mesh->Nverts+1] = V1y;
  EX[1*mesh->Nverts+2] = V4x;  EY[1*mesh->Nverts+2] = V4y;

  EX[2*mesh->Nverts+0] = V3x;  EY[2*mesh->Nverts+0] = V3y;
  EX[2*mesh->Nverts+1] = V2x;  EY[2*mesh->Nverts+1] = V2y;
  EX[2*mesh->Nverts+2] = V5x;  EY[2*mesh->Nverts+2] = V5y;

  EX[3*mesh->Nverts+0] = V1x;  EY[3*mesh->Nverts+0] = V1y;
  EX[3*mesh->Nverts+1] = V3x;  EY[3*mesh->Nverts+1] = V3y;
  EX[3*mesh->Nverts+2] = V6x;  EY[3*mesh->Nverts+2] = V6y;


  dfloat *x = (dfloat *) calloc(NpatchElements*mesh->Np,sizeof(dfloat));
  dfloat *y = (dfloat *) calloc(NpatchElements*mesh->Np,sizeof(dfloat));

  int cnt =0;
  for (int e=0;e<NpatchElements;e++) {
    iint id = e*mesh->Nverts;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    for(iint n=0;n<mesh->Np;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n];
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      x[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      y[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      ++cnt;
    }
  }

  //fake vmap to map facenodes of reference element to face 0 of a neighbour
  int *refFaceMap = (int *) calloc(mesh->Nfaces*mesh->Nfp,sizeof(int));

  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nfp;n++) {
      int id = mesh->faceNodes[f*mesh->Nfp+n];
      dfloat xM = x[id];
      dfloat yM = y[id];

      int matchIndex = mesh->faceNodes[0];
      dfloat xP = x[(f+1)*mesh->Np+matchIndex];
      dfloat yP = y[(f+1)*mesh->Np+matchIndex];
      dfloat mindist = pow(xM-xP,2) + pow(yM-yP,2);
      for (int m=1;m<mesh->Nfp;m++) {
        int idP = mesh->faceNodes[m];
        xP = x[(f+1)*mesh->Np+idP];
        yP = y[(f+1)*mesh->Np+idP];

        dfloat dist = pow(xM-xP,2) + pow(yM-yP,2);
        if(dist<mindist){
          mindist = dist;
          matchIndex = idP;
        }
      }
      if(mindist>1e-3) printf("arggh - bad match, reference patch: x,y=%g,%g\n", xM,yM);
      refFaceMap[f*mesh->Nfp+n] = matchIndex;
    }
  }
  free(EX); free(EY);
  free(x); free(y);



  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Nelements*patchNp*patchNp, sizeof(dfloat));

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){
    dfloat *patchA = patchesInvA[0] + eM*patchNp*patchNp;

    //build the patch A matrix for this element
    BuildFullPatchAx(mesh, basis, tau, lambda, BCType, MS, refFaceMap, eM, patchA);

    // in place inverse (patchA points into patchesInvA[0])
    matrixInverse(patchNp, patchA);
  }

  free(refFaceMap);
  free(MS);
}

void ellipticBuildApproxPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
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

  //We need the halo element's EToB flags to make the patch matrices
  if (mesh->totalHaloPairs) {
    mesh->EToB = (int *) realloc(mesh->EToB,(mesh->Nelements+mesh->totalHaloPairs)*sizeof(int));
    iint *idSendBuffer = (int *) calloc(mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, sizeof(int), mesh->EToB, idSendBuffer, mesh->EToB + mesh->Nelements);
    free(idSendBuffer);
  }

  int NpatchElements = mesh->Nfaces+1;
  int patchNp = mesh->Np*NpatchElements;

 //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x = -1., V4x =  1., V5x = 1., V6x = -3.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y = -3., V5y = 1., V6y =  1.;

  dfloat *EX = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));
  dfloat *EY = (dfloat *) calloc(mesh->Nverts*NpatchElements,sizeof(dfloat));

  EX[0*mesh->Nverts+0] = V1x;  EY[0*mesh->Nverts+0] = V1y;
  EX[0*mesh->Nverts+1] = V2x;  EY[0*mesh->Nverts+1] = V2y;
  EX[0*mesh->Nverts+2] = V3x;  EY[0*mesh->Nverts+2] = V3y;

  EX[1*mesh->Nverts+0] = V2x;  EY[1*mesh->Nverts+0] = V2y;
  EX[1*mesh->Nverts+1] = V1x;  EY[1*mesh->Nverts+1] = V1y;
  EX[1*mesh->Nverts+2] = V4x;  EY[1*mesh->Nverts+2] = V4y;

  EX[2*mesh->Nverts+0] = V3x;  EY[2*mesh->Nverts+0] = V3y;
  EX[2*mesh->Nverts+1] = V2x;  EY[2*mesh->Nverts+1] = V2y;
  EX[2*mesh->Nverts+2] = V5x;  EY[2*mesh->Nverts+2] = V5y;

  EX[3*mesh->Nverts+0] = V1x;  EY[3*mesh->Nverts+0] = V1y;
  EX[3*mesh->Nverts+1] = V3x;  EY[3*mesh->Nverts+1] = V3y;
  EX[3*mesh->Nverts+2] = V6x;  EY[3*mesh->Nverts+2] = V6y;

  dfloat *x = (dfloat *) calloc(NpatchElements*mesh->Np,sizeof(dfloat));
  dfloat *y = (dfloat *) calloc(NpatchElements*mesh->Np,sizeof(dfloat));

  int cnt =0;
  for (int e=0;e<NpatchElements;e++) {
    iint id = e*mesh->Nverts+0;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    for(iint n=0;n<mesh->Np;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[n];
      dfloat sn = mesh->s[n];

      /* physical coordinate of interpolation node */
      x[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      y[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      ++cnt;
    }
  }

  //fake vmap to map facenodes of reference element to face 0 of a neighbour
  int *refFaceMap = (int *) calloc(mesh->Nfaces*mesh->Nfp,sizeof(int));

  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nfp;n++) {
      int id = mesh->faceNodes[f*mesh->Nfp+n];
      dfloat xM = x[id];
      dfloat yM = y[id];

      int matchIndex = mesh->faceNodes[0];
      dfloat xP = x[(f+1)*mesh->Np+matchIndex];
      dfloat yP = y[(f+1)*mesh->Np+matchIndex];
      dfloat mindist = pow(xM-xP,2) + pow(yM-yP,2);
      for (int m=1;m<mesh->Nfp;m++) {
        int idP = mesh->faceNodes[m];
        xP = x[(f+1)*mesh->Np+idP];
        yP = y[(f+1)*mesh->Np+idP];

        dfloat dist = pow(xM-xP,2) + pow(yM-yP,2);
        if(dist<mindist){
          mindist = dist;
          matchIndex = idP;
        }
      }
      if(mindist>1e-3) printf("arggh - bad match, reference patch: x,y=%g,%g\n", xM,yM);
      refFaceMap[f*mesh->Nfp+n] = matchIndex;
    }
  }
  free(EX); free(EY);
  free(x); free(y);



  iint Nperm = pow(mesh->Nfaces,mesh->Nfaces);//all possible configureation of neighbours
  (*Npatches) = Nperm;
  int refPatches = 0;

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(Nperm*patchNp*patchNp, sizeof(dfloat));
  *patchesIndex = (iint*) calloc(mesh->Nelements, sizeof(iint));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(patchNp*patchNp, sizeof(dfloat));


  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildFullPatchAx(mesh, basis, tau, lambda, BCType, MS, refFaceMap, -1, refPatchInvA);
  matrixInverse(patchNp, refPatchInvA);

  //store the permutations of the reference patch
  int blkCounter = 1;
  int *permIndex = (int*) calloc(patchNp, sizeof(int));
  for(iint blk=1;blk<Nperm;++blk){
    iint f0 = blk%mesh->Nfaces;
    iint f1 = (blk/mesh->Nfaces)%mesh->Nfaces;
    iint f2= (blk/(mesh->Nfaces*mesh->Nfaces));

    iint r0 = (f0==0) ? 0: mesh->Nfaces-f0;
    iint r1 = (f1==0) ? 0: mesh->Nfaces-f1;
    iint r2 = (f2==0) ? 0: mesh->Nfaces-f2;

    for(iint n=0;n<mesh->Np;++n){
      permIndex[n+0*mesh->Np] = 0*mesh->Np + n;
      permIndex[n+1*mesh->Np] = 1*mesh->Np + mesh->rmapP[r0*mesh->Np+n];
      permIndex[n+2*mesh->Np] = 2*mesh->Np + mesh->rmapP[r1*mesh->Np+n];
      permIndex[n+3*mesh->Np] = 3*mesh->Np + mesh->rmapP[r2*mesh->Np+n];
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
  for(iint eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildFullPatchAx(mesh, basis, tau, lambda, BCType, MS, refFaceMap, eM, patchA);

    iint eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    iint eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    iint eP2 = mesh->EToE[eM*mesh->Nfaces+2];

    iint fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    iint fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    iint fP2 = mesh->EToF[eM*mesh->Nfaces+2];

    if(eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interiour patch
      iint blk = fP0 + mesh->Nfaces*fP1 + mesh->Nfaces*mesh->Nfaces*fP2;

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

      printf("Element %d's conditioned patch reports cond = %g and rate = %g \n", eM, cond, rate);

      if (rate < 0.7) {
        (*patchesIndex)[eM] = blk;
        refPatches++;
        continue;
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

    (*patchesIndex)[eM] = (*Npatches)-1;
  }

  printf("saving %d full patches\n",*Npatches);
  printf("using %d reference patches\n", refPatches);

  free(patchA); free(invRefAA);
  free(MS); free(refFaceMap);
}
