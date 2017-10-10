#include "ellipticTri2D.h"

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the ipdg patch A matrix for element eM
void BuildLocalIpdgPatchAx(solver_t* solver, mesh2D* mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A);

//returns the BRdg patch A matrix for element eM
void BuildLocalBRdgPatchAx(solver_t* solver, mesh2D* mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A);


void ellipticBuildLocalPatchesTri2D(solver_t* solver, mesh2D* mesh, iint basisNp, dfloat *basis,
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
  dfloat V1x = -1., V2x = 1., V3x =        0.;
  dfloat V1y =  0., V2y = 0., V3y =  sqrt(3.);

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;  refMesh->EY[0] = V1y;
  refMesh->EX[1] = V2x;  refMesh->EY[1] = V2y;
  refMesh->EX[2] = V3x;  refMesh->EY[2] = V3y;

  refMesh->EToV = (iint*) calloc(mesh->Nverts, sizeof(iint));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;

  refMesh->EToB = (iint*) calloc(mesh->Nfaces,sizeof(iint));
  for (iint n=0;n<mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTri2D(refMesh, mesh->N);
  meshPhysicalNodesTri2D(refMesh);
  meshGeometricFactorsTri2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsTri2D(refMesh);

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  if (strstr(options,"IPDG")) {
    BuildLocalIpdgPatchAx(solver, refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  } else if (strstr(options,"BRDG")) {
    BuildLocalBRdgPatchAx(solver, refMesh, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  }
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(iint eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    if (strstr(options,"IPDG")) {
      BuildLocalIpdgPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, eM, patchA);
    } else if (strstr(options,"BRDG")) {
      BuildLocalBRdgPatchAx(solver, mesh, basis, tau, lambda, BCType, MS, eM, patchA);
    }


    iint eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    iint eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    iint eP2 = mesh->EToE[eM*mesh->Nfaces+2];

    iint fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    iint fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    iint fP2 = mesh->EToF[eM*mesh->Nfaces+2];

    if(eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interior patch

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

  free(refMesh);
  free(patchA); free(invRefAA);
  free(MS);
}


//returns the ipdg patch A matrix for element eM
void BuildLocalIpdgPatchAx(solver_t* solver, mesh2D* mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A) {

  iint vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];
  dfloat J = mesh->vgeo[vbase+JID];

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
    iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
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

//returns the ipdg patch A matrix for element eM
void BuildLocalBRdgPatchAx(solver_t* solver, mesh2D* mesh, dfloat *basis, dfloat tau, dfloat lambda, iint* BCType,
                        dfloat *MS, iint eM, dfloat *A) {

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;

  iint vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];
  dfloat J = mesh->vgeo[vbase+JID];

  dfloat* Gx = (dfloat*) calloc(mesh->Np*mesh->Np*(Nfaces+1),sizeof(dfloat));
  dfloat* Gy = (dfloat*) calloc(mesh->Np*mesh->Np*(Nfaces+1),sizeof(dfloat));

  dfloat *qmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *qmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *QmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *QmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));

  dfloat* Ae = (dfloat*) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  for(iint n=0;n<Np;++n){
    for(iint m=0;m<Np;++m){
      Gx[m+n*Np] = drdx*mesh->Dr[m+n*Np]+dsdx*mesh->Ds[m+n*Np];
      Gy[m+n*Np] = drdy*mesh->Dr[m+n*Np]+dsdy*mesh->Ds[m+n*Np];
    }
  }

  for (iint fM=0;fM<Nfaces;fM++) {
    // load surface geofactors for this face
    iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat invJ = mesh->sgeo[sid+IJID];

    iint eP = mesh->EToE[eM*Nfaces+fM];
    
    for (iint m=0;m<Np;m++) {
      // extract trace nodes
      for (iint i=0;i<Nfp;i++) {
        // double check vol geometric factors are in halo storage of vgeo
        iint idM    = eM*Nfp*Nfaces+fM*Nfp+i;
        iint vidM   = mesh->faceNodes[i+fM*Nfp];

        qmM[i] =0;
        if (vidM == m) qmM[i] =1;
      }

      int bcD = 0;
      int bcN = 0;
      if (eP < 0) {
        int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
        iint bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
        if(bcType==1){ // Dirichlet
          bcD = 1;
          bcN = 0;
        } else if (bcType==2){ // Neumann
          bcD = 0;
          bcN = 1;
        } else { // Neumann for now
          bcD = 0;
          bcN = 1;
        }
      }

      for (iint n=0;n<Np;n++) {
        for (iint i=0;i<Nfp;i++) {
          Gx[m+n*Np] += -0.5*(1-bcN)*(1+bcD)*sJ*invJ*nx*mesh->LIFT[i+fM*Nfp+n*Nfp*Nfaces]*qmM[i];
          Gy[m+n*Np] += -0.5*(1-bcN)*(1+bcD)*sJ*invJ*ny*mesh->LIFT[i+fM*Nfp+n*Nfp*Nfaces]*qmM[i];
        }
      }
    }

    if (eP>-1) {
      int fP = mesh->EToF[eM*Nfaces+fM];
      iint sidP = mesh->Nsgeo*(eP*Nfaces+fP);
      dfloat nxP = mesh->sgeo[sidP+NXID];
      dfloat nyP = mesh->sgeo[sidP+NYID];
      dfloat sJP = mesh->sgeo[sidP+SJID];
      dfloat invJP = mesh->sgeo[sidP+IJID];

      for (iint m=0;m<Np;m++) {
        // extract trace nodes
        for (iint i=0;i<Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          iint idM    = eP*Nfp*Nfaces+fP*Nfp+i;
          iint vidP   = mesh->vmapP[idM]%Np; // only use this to identify location of positive trace vgeo

          qmP[i] =0;
          if (vidP == m) qmP[i] =1;
        }

        for (iint n=0;n<Np;n++) {
          for (iint i=0;i<Nfp;i++) {
            Gx[m+n*Np+(fM+1)*Np*Np] += 0.5*sJP*invJP*nxP*mesh->LIFT[i+fP*Nfp+n*Nfp*Nfaces]*qmP[i];
            Gy[m+n*Np+(fM+1)*Np*Np] += 0.5*sJP*invJP*nyP*mesh->LIFT[i+fP*Nfp+n*Nfp*Nfaces]*qmP[i];
          }
        }
      }
    }
  }

  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      for (int k=0;k<Np;k++) {
        Ae[m+n*Np] += (drdx*mesh->Dr[k+n*Np]+dsdx*mesh->Ds[k+n*Np])*Gx[m+k*Np];
        Ae[m+n*Np] += (drdy*mesh->Dr[k+n*Np]+dsdy*mesh->Ds[k+n*Np])*Gy[m+k*Np];
      }  
    }
    Ae[n+n*Np] -= lambda;  
  }

  for (iint m=0;m<Np;m++) {
    for (iint fM=0;fM<Nfaces;fM++) {
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat invJ = mesh->sgeo[sid+IJID];

      iint eP = mesh->EToE[eM*Nfaces+fM];
      if (eP<0) eP = eM;
    
      // extract trace matrix from Gx and Gy
      for (iint i=0;i<Nfp;i++) {
        // double check vol geometric factors are in halo storage of vgeo
        iint idM    = eM*Nfp*Nfaces+fM*Nfp+i;
        iint vidM   = mesh->faceNodes[i+fM*Nfp];
        iint vidP   = mesh->vmapP[idM]%Np; // only use this to identify location of positive trace vgeo

        qmM[i] = 0;          
        if (vidM == m) qmM[i] = 1;
        QmM[i] = nx*Gx[m+vidM*Np]+ny*Gy[m+vidM*Np];
        QmP[i] = nx*Gx[m+vidP*Np+(fM+1)*Np*Np] + ny*Gy[m+vidP*Np+(fM+1)*Np*Np];                
      }

      int bcD = 0;
      int bcN = 0;
      eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) {
        int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
        iint bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
        if(bcType==1){ // Dirichlet
          bcD = 1;
          bcN = 0;
        } else if (bcType==2){ // Neumann
          bcD = 0;
          bcN = 1;
        } else { // Neumann for now
          bcD = 0;
          bcN = 1;
        }
      }

      for (int n=0;n<Np;n++) {
        for (iint i=0;i<Nfp;i++) {
          Ae[m+n*Np] += -0.5*(1+bcN)*(1-bcD)*sJ*invJ*mesh->LIFT[i+fM*Nfp+n*Nfp*Nfaces]*QmM[i];
          Ae[m+n*Np] +=  0.5*(1-bcN)*(1-bcD)*sJ*invJ*mesh->LIFT[i+fM*Nfp+n*Nfp*Nfaces]*QmP[i];
          Ae[m+n*Np] += -0.5*(1-bcN)*(1+bcD)*sJ*invJ*mesh->LIFT[i+fM*Nfp+n*Nfp*Nfaces]*tau*qmM[i];
        }
      }
    }
  }

  //multiply by mass matrix 
  for (int n=0;n<Np;n++) {
    for (int m=0;m<Np;m++) {
      dfloat Anm = 0.;
      for (int k=0;k<Np;k++) {
        Anm += mesh->MM[k+n*Np]*Ae[m+k*Np]; 
      }
      Anm *= J;

      A[m+n*Np] = -Anm;
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

  free(Gx); free(Gy);
  free(qmM); free(qmP);
  free(QmM); free(QmP);
  free(Ae);
}

