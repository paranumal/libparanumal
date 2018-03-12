#include "ellipticHex3D.h"

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the patch A matrix for element eM
void BuildLocalIpdgPatchAx(mesh3D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *A);

//returns the C0FEM patch A matrix for element eM
void BuildLocalContinuousPatchAx(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dfloat *A);

void ellipticBuildLocalPatchesHex3D(solver_t *solver, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   dlong *Npatches, dlong **patchesIndex, dfloat **patchesInvA,
                                   const char *options){

  mesh3D *mesh = solver->mesh;

  // build some monolithic basis arrays
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nk=0;nk<mesh->N+1;++nk){
    for(int nj=0;nj<mesh->N+1;++nj){
      for(int ni=0;ni<mesh->N+1;++ni){

        int node = 0;

        for(int k=0;k<mesh->N+1;++k){
          for(int j=0;j<mesh->N+1;++j){
            for(int i=0;i<mesh->N+1;++i){

              if(nk==k && nj==j && ni==i)
                B[mode*mesh->Np+node] = 1;
              if(nj==j && nk==k)
                Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
              if(ni==i && nk==k)
                Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
              if(ni==i && nj==j)
                Bt[mode*mesh->Np+node] = mesh->D[nk+mesh->Nq*k]; 
              
              ++node;
            }
          }
        }
        
        ++mode;
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
  mesh3D *refMesh = (mesh3D*) calloc(1,sizeof(mesh3D));
  memcpy(refMesh,mesh,sizeof(mesh3D));

  //vertices of reference patch
  dfloat V1x = -1., V2x =  1., V3x =  1., V4x = -1., V5x = -1., V6x =  1., V7x =  1., V8x = -1.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y =  1., V5y = -1., V6y = -1., V7y =  1., V8y =  1.;
  dfloat V1z = -1., V2z = -1., V3z = -1., V4z = -1., V5z =  1., V6z =  1., V7z =  1., V8z =  1.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EZ = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;  refMesh->EY[0] = V1y; refMesh->EZ[0] = V1z;
  refMesh->EX[1] = V2x;  refMesh->EY[1] = V2y; refMesh->EZ[1] = V2z;
  refMesh->EX[2] = V3x;  refMesh->EY[2] = V3y; refMesh->EZ[2] = V3z;
  refMesh->EX[3] = V4x;  refMesh->EY[3] = V4y; refMesh->EZ[3] = V4z;
  refMesh->EX[4] = V5x;  refMesh->EY[4] = V5y; refMesh->EZ[4] = V5z;
  refMesh->EX[5] = V6x;  refMesh->EY[5] = V6y; refMesh->EZ[5] = V6z;
  refMesh->EX[6] = V7x;  refMesh->EY[6] = V7y; refMesh->EZ[6] = V7z;
  refMesh->EX[7] = V8x;  refMesh->EY[7] = V8y; refMesh->EZ[7] = V8z;

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;
  refMesh->EToV[4] = 4;
  refMesh->EToV[5] = 5;
  refMesh->EToV[6] = 6;
  refMesh->EToV[7] = 7;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n=0;n<mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesHex3D(refMesh, mesh->N);
  meshPhysicalNodesHex3D(refMesh);
  meshGeometricFactorsHex3D(refMesh);
  meshConnectFaceNodes3D(refMesh);
  meshSurfaceGeometricFactorsHex3D(refMesh);

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildLocalIpdgPatchAx(refMesh, tau, lambda, BCType, B,Br,Bs,Bt, 0, refPatchInvA);
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildLocalIpdgPatchAx(mesh, tau, lambda, BCType, B,Br,Bs,Bt, eM, patchA);

    dlong eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    dlong eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    dlong eP2 = mesh->EToE[eM*mesh->Nfaces+2];
    dlong eP3 = mesh->EToE[eM*mesh->Nfaces+3];
    dlong eP4 = mesh->EToE[eM*mesh->Nfaces+4];
    dlong eP5 = mesh->EToE[eM*mesh->Nfaces+5];

    int fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    int fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    int fP2 = mesh->EToF[eM*mesh->Nfaces+2];
    int fP3 = mesh->EToF[eM*mesh->Nfaces+3];
    int fP4 = mesh->EToF[eM*mesh->Nfaces+4];
    int fP5 = mesh->EToF[eM*mesh->Nfaces+5];

    if(eP0>=0 && eP1>=0 && eP2>=0 && eP3>=0 && eP4>=0 && eP5>=0){ //check if this is an interior patch

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

  printf("saving "dlongFormat" full patches\n",*Npatches);
  printf("using "dlongFormat" reference patches\n", refPatches);

  free(refMesh);
  free(patchA); free(invRefAA);
  free(B); free(Br); free(Bs);
}


//returns the patch A matrix for element eM
void BuildLocalIpdgPatchAx(mesh3D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *A) {

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      A[n*mesh->Np+m] = 0;

      // (grad phi_n, grad phi_m)_{D^e}
      for(int i=0;i<mesh->Np;++i){
        dlong base = eM*mesh->Np*mesh->Nvgeo + i;
        dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
        dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
        dfloat drdz = mesh->vgeo[base+mesh->Np*RZID];
        dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
        dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
        dfloat dsdz = mesh->vgeo[base+mesh->Np*SZID];
        dfloat dtdx = mesh->vgeo[base+mesh->Np*TXID];
        dfloat dtdy = mesh->vgeo[base+mesh->Np*TYID];
        dfloat dtdz = mesh->vgeo[base+mesh->Np*TZID];
        dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
        
        int idn = n*mesh->Np+i;
        int idm = m*mesh->Np+i;
        dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
        dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
        dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];    
        dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm] + dtdx*Bt[idm];
        dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm] + dtdy*Bt[idm];
        dfloat dlmdz = drdz*Br[idm] + dsdz*Bs[idm] + dtdz*Bt[idm];
        A[n*mesh->Np+m] += JW*(dlndx*dlmdx+dlndy*dlmdy+dlndz*dlmdz);
        A[n*mesh->Np+m] += lambda*JW*B[idn]*B[idm];
      }

      for (int fM=0;fM<mesh->Nfaces;fM++) {
        // accumulate flux terms for negative and positive traces
        for(int i=0;i<mesh->Nfp;++i){
          int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

          // grab vol geofacs at surface nodes
          dlong baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
          dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
          dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
          dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
          dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
          dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
          dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];
          dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
          dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
          dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];

          // grab surface geometric factors
          dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
          dfloat nx = mesh->sgeo[base+NXID];
          dfloat ny = mesh->sgeo[base+NYID];
          dfloat nz = mesh->sgeo[base+NZID];
          dfloat wsJ = mesh->sgeo[base+WSJID];
          dfloat hinv = mesh->sgeo[base+IHID];
          
          // form negative trace terms in IPDG
          int idnM = n*mesh->Np+vidM; 
          int idmM = m*mesh->Np+vidM;

          dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
          dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
          dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
          dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
          dfloat lnM = B[idnM];

          dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM] + dtdxM*Bt[idmM];
          dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM] + dtdyM*Bt[idmM];
          dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM] + dtdzM*Bt[idmM];
          dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM+nz*dlmdzM;
          dfloat lmM = B[idmM];
          
          dfloat penalty = tau*hinv;     
          int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag

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

          A[n*mesh->Np+m] += -0.5*(1+bcD)*(1-bcN)*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
          A[n*mesh->Np+m] += -0.5*(1+bcD)*(1-bcN)*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
          A[n*mesh->Np+m] += +0.5*(1+bcD)*(1-bcN)*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)
        }
      }
    }
  }
}

void BuildLocalContinuousPatchAx(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat* Bs, dfloat* Bt, dfloat *A) {

  mesh3D *mesh = solver->mesh;

  for (int nz=0;nz<mesh->Nq;nz++) {
  for (int ny=0;ny<mesh->Nq;ny++) {
  for (int nx=0;nx<mesh->Nq;nx++) {
    int idn = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
    if (solver->mapB[idn+eM*mesh->Np]!=1) {
      for (int mz=0;mz<mesh->Nq;mz++) {
      for (int my=0;my<mesh->Nq;my++) {
      for (int mx=0;mx<mesh->Nq;mx++) {  
        int idm = mx+my*mesh->Nq+mz*mesh->Nq*mesh->Nq;
        int iid = idn*mesh->Np + idm;
        if (solver->mapB[idm+eM*mesh->Np]==1) continue;
            
        int id;
        A[iid] = 0;

        if ((ny==my)&&(nz==mz)) {
          for (int k=0;k<mesh->Nq;k++) {
            id = k+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
            dfloat Grr = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];

            A[iid] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[mx+k*mesh->Nq];
          }
        }

        if (nz==mz) {
          id = mx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
          dfloat Grs = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
          A[iid] += Grs*mesh->D[nx+mx*mesh->Nq]*mesh->D[my+ny*mesh->Nq];

          id = nx+my*mesh->Nq+nz*mesh->Nq*mesh->Nq;
          dfloat Gsr = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
          A[iid] += Gsr*mesh->D[mx+nx*mesh->Nq]*mesh->D[ny+my*mesh->Nq];
        }

        if (ny==my) {
          id = mx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
          dfloat Grt = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
          A[iid] += Grt*mesh->D[nx+mx*mesh->Nq]*mesh->D[mz+nz*mesh->Nq];

          id = nx+ny*mesh->Nq+mz*mesh->Nq*mesh->Nq;
          dfloat Gst = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G02ID*mesh->Np];
          A[iid] += Gst*mesh->D[mx+nx*mesh->Nq]*mesh->D[nz+mz*mesh->Nq];
        }

        if ((nx==mx)&&(nz==mz)) {
          for (int k=0;k<mesh->Nq;k++) {
            id = nx+k*mesh->Nq+nz*mesh->Nq*mesh->Nq;
            dfloat Gss = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];

            A[iid] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[my+k*mesh->Nq];
          }
        }
        
        if (nx==mx) {
          id = nx+my*mesh->Nq+nz*mesh->Nq*mesh->Nq;
          dfloat Gst = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
          A[iid] += Gst*mesh->D[ny+my*mesh->Nq]*mesh->D[mz+nz*mesh->Nq];

          id = nx+ny*mesh->Nq+mz*mesh->Nq*mesh->Nq;
          dfloat Gts = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G12ID*mesh->Np];
          A[iid] += Gts*mesh->D[my+ny*mesh->Nq]*mesh->D[nz+mz*mesh->Nq];
        }

        if ((nx==mx)&&(ny==my)) {
          for (int k=0;k<mesh->Nq;k++) {
            id = nx+ny*mesh->Nq+k*mesh->Nq*mesh->Nq;
            dfloat Gtt = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G22ID*mesh->Np];

            A[iid] += Gtt*mesh->D[nz+k*mesh->Nq]*mesh->D[mz+k*mesh->Nq];
          }
        }
        
        if ((nx==mx)&&(ny==my)&&(nz==mz)) {
          id = nx + ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
          dfloat JW = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + GWJID*mesh->Np];
          A[iid] += JW*lambda;
        }
      }
      }
      }
    } else {
      int iid = idn*mesh->Np + idn;
      A[iid] = 1; //just put a 1 so A is invertable
    }
  }
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
