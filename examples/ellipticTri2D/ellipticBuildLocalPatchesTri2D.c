#include "ellipticTri2D.h"

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the ipdg patch A matrix for element eM
void BuildLocalIpdgPatchAx(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A);

//returns the BRdg patch A matrix for element eM
void BuildLocalBRdgPatchAx(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A);

//returns the C0FEM patch A matrix for element eM
void BuildLocalContinuousPatchAx(solver_t* solver, mesh2D* mesh, dfloat lambda,
                                  dlong eM, dfloat *Srr, dfloat *Srs, 
                                  dfloat *Sss, dfloat *MM, dfloat *A);

void ellipticBuildLocalPatchesTri2D(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis,
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

  //patch inverse storage
  *patchesInvA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  *patchesIndex = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  //temp patch storage
  dfloat *patchA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *invRefAA = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  (*Npatches) = 1;
  dlong refPatches = 0;


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

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n=0;n<mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesTri2D(refMesh, mesh->N);
  meshPhysicalNodesTri2D(refMesh);
  meshGeometricFactorsTri2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsTri2D(refMesh);

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  if (strstr(options,"IPDG")) {
    BuildLocalIpdgPatchAx(solver, refMesh, basisNp, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  } else if (strstr(options,"BRDG")) {
    BuildLocalBRdgPatchAx(solver, refMesh, basisNp, basis, tau, lambda, BCType, MS, 0, refPatchInvA);
  }
#if 0
  for (int n=0;n<mesh->Np;n++) {
    for (int m=0;m<mesh->Np;m++) {
      printf("%4.2f \t", refPatchInvA[m+n*mesh->Np]);
    }
    printf("\n");
  }
  printf("\n");
#endif
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    if (strstr(options,"IPDG")) {
      BuildLocalIpdgPatchAx(solver, mesh, basisNp, basis, tau, lambda, BCType, MS, eM, patchA);
    } else if (strstr(options,"BRDG")) {
      BuildLocalBRdgPatchAx(solver, mesh, basisNp, basis, tau, lambda, BCType, MS, eM, patchA);
    }
#if 0
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
        printf("%4.2f \t", patchA[m+n*mesh->Np]);
      }
      printf("\n");
    }
    printf("\n");
#endif
    dlong eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    dlong eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    dlong eP2 = mesh->EToE[eM*mesh->Nfaces+2];

    int fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    int fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    int fP2 = mesh->EToF[eM*mesh->Nfaces+2];

    if(eP0>=0 && eP1>=0 && eP2>=0){ //check if this is an interior patch

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
        dlong id = ((*Npatches)-1)*mesh->Np*mesh->Np + n*mesh->Np + m;
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
void BuildLocalIpdgPatchAx(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A) {

  dlong vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];
  dfloat J = mesh->vgeo[vbase+JID];

  dfloat *Ae = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      Ae[n*mesh->Np+m]  = J*lambda*mesh->MM[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*drdx*drdx*mesh->Srr[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*drdx*dsdx*mesh->Srs[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*dsdx*drdx*mesh->Ssr[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*dsdx*dsdx*mesh->Sss[n*mesh->Np+m];

      Ae[n*mesh->Np+m] += J*drdy*drdy*mesh->Srr[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*drdy*dsdy*mesh->Srs[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*dsdy*drdy*mesh->Ssr[n*mesh->Np+m];
      Ae[n*mesh->Np+m] += J*dsdy*dsdy*mesh->Sss[n*mesh->Np+m];
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (solver->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        Ae[n*mesh->Np+m] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
      }
    }
  }

  for (int fM=0;fM<mesh->Nfaces;fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat hinv = mesh->sgeo[sid+IHID];

    int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

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

        // OP11 = OP11 + 0.5*( gtau*mmE )
        dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
        Ae[nM*mesh->Np+mM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
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

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          Ae[nM*mesh->Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          Ae[nM*mesh->Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
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

          // OP11 = OP11 + (- Dn1'*mmE );
          Ae[n*mesh->Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
          Ae[n*mesh->Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
        }
      }
    }
  }

  for(int j=0;j<basisNp;++j){
    for(int i=0;i<basisNp;++i){
      dfloat val = 0;
      for (int n=0;n<mesh->Np;n++) {
        for (int m=0;m<mesh->Np;m++) {
          val += basis[n*basisNp+j]*Ae[n*mesh->Np+m]*basis[m*basisNp+i];
        }
      }

      A[i+j*basisNp] = val;
    }
  }

  free(Ae);
}

//returns the ipdg patch A matrix for element eM
void BuildLocalBRdgPatchAx(solver_t* solver, mesh2D* mesh, int basisNp, dfloat *basis, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *MS, dlong eM, dfloat *A) {

  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;

  /* Construct gradient as block matrix and load it to the halo */
  dfloat  *Gx = (dfloat *) calloc(Np*Np*(Nfaces+1),sizeof(dfloat));
  dfloat  *Gy = (dfloat *) calloc(Np*Np*(Nfaces+1),sizeof(dfloat));

  dlong vbase = eM*mesh->Nvgeo;
  dfloat drdx = mesh->vgeo[vbase+RXID];
  dfloat drdy = mesh->vgeo[vbase+RYID];
  dfloat dsdx = mesh->vgeo[vbase+SXID];
  dfloat dsdy = mesh->vgeo[vbase+SYID];
  dfloat J = mesh->vgeo[vbase+JID];

  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      Gx[m+n*Np] = drdx*mesh->Dr[m+n*Np]+dsdx*mesh->Ds[m+n*Np];
      Gy[m+n*Np] = drdy*mesh->Dr[m+n*Np]+dsdy*mesh->Ds[m+n*Np];
    }
  }

  for (int fM=0;fM<Nfaces;fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat invJ = mesh->sgeo[sid+IJID];

    dlong eP = mesh->EToE[eM*Nfaces+fM];
    int fP = mesh->EToF[eM*Nfaces+fM];
    dfloat sw = 1.f; //guard against unconnected elements (happens in reference patch)
    if (eP < 0) {eP = eM; sw = 0;}
    if (fP < 0) fP = fM;

    // load surface geofactors for neighbor's face
    dlong sidP = mesh->Nsgeo*(eP*Nfaces+fP);
    dfloat nxP = mesh->sgeo[sidP+NXID];
    dfloat nyP = mesh->sgeo[sidP+NYID];
    dfloat sJP = mesh->sgeo[sidP+SJID];
    dfloat invJP = mesh->sgeo[sidP+IJID];

    int bcD = 0, bcN =0;
    int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
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

        dlong idM = eP*Nfp*Nfaces+fP*Nfp+m;
        int mP = (int) (mesh->vmapP[idM]%Np);

        dfloat LIFTfnmM = sJ*invJ*mesh->LIFT[m + fM*Nfp + n*Nfp*Nfaces];
        dfloat LIFTfnmP = sJP*invJP*mesh->LIFT[m + fP*Nfp + n*Nfp*Nfaces];

        // G = sJ/J*LIFT*n*[[ uP-uM ]]
        Gx[mM+n*Np] += -0.5*(1-bcN)*(1+bcD)*nx*LIFTfnmM;
        Gy[mM+n*Np] += -0.5*(1-bcN)*(1+bcD)*ny*LIFTfnmM;

        Gx[mP+n*Np+(fM+1)*Np*Np] +=  0.5*sw*(1-bcN)*(1-bcD)*nxP*LIFTfnmP;
        Gy[mP+n*Np+(fM+1)*Np*Np] +=  0.5*sw*(1-bcN)*(1-bcD)*nyP*LIFTfnmP;
      }
    }
  }

  dfloat *Ae = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  /* start with stiffness matrix  */
  for(int n=0;n<Np;++n){
    for(int m=0;m<Np;++m){
      Ae[n*Np+m]  = J*lambda*mesh->MM[n*Np+m];
      Ae[n*Np+m] += J*drdx*drdx*mesh->Srr[n*Np+m];
      Ae[n*Np+m] += J*drdx*dsdx*mesh->Srs[n*Np+m];
      Ae[n*Np+m] += J*dsdx*drdx*mesh->Ssr[n*Np+m];
      Ae[n*Np+m] += J*dsdx*dsdx*mesh->Sss[n*Np+m];

      Ae[n*Np+m] += J*drdy*drdy*mesh->Srr[n*Np+m];
      Ae[n*Np+m] += J*drdy*dsdy*mesh->Srs[n*Np+m];
      Ae[n*Np+m] += J*dsdy*drdy*mesh->Ssr[n*Np+m];
      Ae[n*Np+m] += J*dsdy*dsdy*mesh->Sss[n*Np+m];
    }
  }

  for (int fM=0;fM<Nfaces;fM++) {
    // load surface geofactors for this face
    dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
    dfloat nx = mesh->sgeo[sid+NXID];
    dfloat ny = mesh->sgeo[sid+NYID];
    dfloat sJ = mesh->sgeo[sid+SJID];
    dfloat hinv = mesh->sgeo[sid+IHID];

    int bcD = 0, bcN =0;
    int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
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
    dfloat *MSf = MS + fM*Nfp*Nfp;

    // penalty term just involves face nodes
    for(int n=0;n<Nfp;++n){
      for(int m=0;m<Nfp;++m){
        int nM = mesh->faceNodes[fM*Nfp+n];
        int mM = mesh->faceNodes[fM*Nfp+m];

        dfloat MSfnm = sJ*MSf[n*Nfp+m];
        Ae[nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*tau*MSfnm;
      }
    }

    // now add differential surface terms
    for(int n=0;n<Nfp;++n){
      for(int m=0;m<Np;++m){
        int nM = mesh->faceNodes[fM*Nfp+n];

        for(int i=0;i<Nfp;++i){
          int iM = mesh->faceNodes[fM*Nfp+i];
          int iP = (int) (mesh->vmapP[i+fM*Nfp+eM*Nfp*Nfaces]%Np);

          dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

          dfloat DxMim = Gx[m+iM*Np];
          dfloat DyMim = Gy[m+iM*Np];

          dfloat DxPim = Gx[m+iP*Np+(fM+1)*Np*Np];
          dfloat DyPim = Gy[m+iP*Np+(fM+1)*Np*Np];

          Ae[m+nM*Np] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          Ae[m+nM*Np] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

          Ae[m+nM*Np] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
          Ae[m+nM*Np] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
        }
      }
    }

    for(int n=0;n<Np;++n){
      for(int m=0;m<Nfp;++m){
        int mM = mesh->faceNodes[fM*Nfp+m];
        int mP = (int) (mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

        for(int i=0;i<Nfp;++i){
          int iM = mesh->faceNodes[fM*Nfp+i];

          dfloat MSfim = sJ*MSf[i*Nfp+m];

          dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
          dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

          Ae[mM+n*Np] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
          Ae[mM+n*Np] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
        }
      }
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (solver->allNeumann) {
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        Ae[n*Np+m] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
      }
    }
  }

  for(int j=0;j<basisNp;++j){
    for(int i=0;i<basisNp;++i){
      dfloat val = 0;
      for (int n=0;n<mesh->Np;n++) {
        for (int m=0;m<mesh->Np;m++) {
          val += basis[n*basisNp+j]*Ae[n*mesh->Np+m]*basis[m*basisNp+i];
        }
      }

      A[i+j*basisNp] = val;
    }
  }

  free(Ae);

  free(Gx); free(Gy);
}

//returns the continuous C0 patch A matrix for element eM
void BuildLocalContinuousPatchAx(solver_t* solver, mesh2D* mesh, dfloat lambda,
                                  dlong eM, dfloat *Srr, dfloat *Srs, 
                                  dfloat *Sss, dfloat *MM, dfloat *A) {

  dlong gbase = eM*mesh->Nggeo;
  dfloat Grr = mesh->ggeo[gbase + G00ID];
  dfloat Grs = mesh->ggeo[gbase + G01ID];
  dfloat Gss = mesh->ggeo[gbase + G11ID];
  dfloat J   = mesh->ggeo[gbase + GWJID];

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    if (solver->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
      for(int m=0;m<mesh->Np;++m){
        if (solver->mapB[m+eM*mesh->Np]!=1) {//dont fill rows for masked nodes
          A[n*mesh->Np+m] = J*lambda*MM[m+n*mesh->Np];
          A[n*mesh->Np+m] += Grr*Srr[m+n*mesh->Np];
          A[n*mesh->Np+m] += Grs*Srs[m+n*mesh->Np];
          A[n*mesh->Np+m] += Gss*Sss[m+n*mesh->Np];
        } else {
          A[n*mesh->Np+m] = 0;
        }
      }
    } else {
      A[n+n*mesh->Np] = 1; //just put a 1 so A is invertable
    }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (solver->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (solver->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        for(int m=0;m<mesh->Np;++m){
          if (solver->mapB[m+eM*mesh->Np]==1) continue; //skip masked nodes
          A[n*mesh->Np+m] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
        }
      } 
    }
  }
}