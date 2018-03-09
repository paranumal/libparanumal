#include "ellipticQuad2D.h"

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

//returns the patch A matrix for element eM
void BuildLocalIpdgPatchAx(mesh2D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, dlong eM, dfloat *A);

//returns the C0FEM patch A matrix for element eM
void BuildLocalContinuousPatchAx(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *A);

void ellipticBuildLocalPatchesQuad2D(solver_t *solver, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   dlong *Npatches, dlong **patchesIndex, dfloat **patchesInvA,
                                   const char *options){

  mesh2D *mesh = solver->mesh;

  // build some monolithic basis arrays
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nj=0;nj<mesh->N+1;++nj){
    for(int ni=0;ni<mesh->N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh->N+1;++j){
        for(int i=0;i<mesh->N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh->Np+node] = 1;
          if(nj==j)
            Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
          if(ni==i)
            Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
          
          ++node;
        }
      }
      ++mode;
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
  dfloat V1x = -1., V2x =  1., V3x =  1., V4x = -1.;
  dfloat V1y = -1., V2y = -1., V3y =  1., V4y =  1.;

  refMesh->Nelements = 1;

  refMesh->EX = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));
  refMesh->EY = (dfloat *) calloc(mesh->Nverts,sizeof(dfloat));

  refMesh->EX[0] = V1x;  refMesh->EY[0] = V1y;
  refMesh->EX[1] = V2x;  refMesh->EY[1] = V2y;
  refMesh->EX[2] = V3x;  refMesh->EY[2] = V3y;
  refMesh->EX[3] = V4x;  refMesh->EY[3] = V4y;

  refMesh->EToV = (hlong*) calloc(mesh->Nverts, sizeof(hlong));

  refMesh->EToV[0] = 0;
  refMesh->EToV[1] = 1;
  refMesh->EToV[2] = 2;
  refMesh->EToV[3] = 3;

  refMesh->EToB = (int*) calloc(mesh->Nfaces,sizeof(int));
  for (int n=0;n<mesh->Nfaces;n++) refMesh->EToB[n] = 0;

  meshConnect(refMesh);
  meshLoadReferenceNodesQuad2D(refMesh, mesh->N);
  meshPhysicalNodesQuad2D(refMesh);
  meshGeometricFactorsQuad2D(refMesh);
  meshConnectFaceNodes2D(refMesh);
  meshSurfaceGeometricFactorsQuad2D(refMesh);

  //start with reference patch
  dfloat *refPatchInvA = *patchesInvA;
  BuildLocalIpdgPatchAx(refMesh, tau, lambda, BCType, B,Br,Bs, 0, refPatchInvA);
  matrixInverse(mesh->Np, refPatchInvA);

  // loop over all elements
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    //build the patch A matrix for this element
    BuildLocalIpdgPatchAx(mesh, tau, lambda, BCType, B,Br,Bs, eM, patchA);

    dlong eP0 = mesh->EToE[eM*mesh->Nfaces+0];
    dlong eP1 = mesh->EToE[eM*mesh->Nfaces+1];
    dlong eP2 = mesh->EToE[eM*mesh->Nfaces+2];
    dlong eP3 = mesh->EToE[eM*mesh->Nfaces+3];

    int fP0 = mesh->EToF[eM*mesh->Nfaces+0];
    int fP1 = mesh->EToF[eM*mesh->Nfaces+1];
    int fP2 = mesh->EToF[eM*mesh->Nfaces+2];
    int fP3 = mesh->EToF[eM*mesh->Nfaces+3];

    if(eP0>=0 && eP1>=0 && eP2>=0 && eP3>=0){ //check if this is an interior patch

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
void BuildLocalIpdgPatchAx(mesh2D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, dlong eM, dfloat *A) {

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    for(int m=0;m<mesh->Np;++m){
      A[n*mesh->Np+m] = 0;

      // (grad phi_n, grad phi_m)_{D^e}
      for(int i=0;i<mesh->Np;++i){
        dlong base = eM*mesh->Np*mesh->Nvgeo + i;
        dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
        dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
        dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
        dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
        dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
        
        int idn = n*mesh->Np+i;
        int idm = m*mesh->Np+i;
        dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
        dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
        dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm];
        dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm];
        A[n*mesh->Np+m] += JW*(dlndx*dlmdx+dlndy*dlmdy);
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
          dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
          dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];

          // grab surface geometric factors
          dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
          dfloat nx = mesh->sgeo[base+NXID];
          dfloat ny = mesh->sgeo[base+NYID];
          dfloat wsJ = mesh->sgeo[base+WSJID];
          dfloat hinv = mesh->sgeo[base+IHID];
          
          // form negative trace terms in IPDG
          int idnM = n*mesh->Np+vidM; 
          int idmM = m*mesh->Np+vidM;

          dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
          dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
          dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
          dfloat lnM = B[idnM];

          dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];
          dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];
          dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM;
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
                                  dlong eM, dfloat *B, dfloat *Br, dfloat* Bs, dfloat *A) {

  mesh2D *mesh = solver->mesh;

  for (int ny=0;ny<mesh->Nq;ny++) {
    for (int nx=0;nx<mesh->Nq;nx++) {
      if (solver->mapB[nx+ny*mesh->Nq+eM*mesh->Np]!=1) continue; //skip masked nodes
      for (int my=0;my<mesh->Nq;my++) {
        for (int mx=0;mx<mesh->Nq;mx++) {
          if (solver->mapB[mx+my*mesh->Nq+eM*mesh->Np]!=1) continue; //skip masked nodes
          
          int id;
          int iid = (nx+ny*mesh->Nq)*mesh->Np + mx+my*mesh->Nq;
          A[iid] = 0;

          if (ny==my) {
            for (int k=0;k<mesh->Nq;k++) {
              id = k+ny*mesh->Nq;
              dfloat Grr = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G00ID*mesh->Np];

              A[iid] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[mx+k*mesh->Nq];
            }
          }

          id = mx+ny*mesh->Nq;
          dfloat Grs = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
          A[iid] += Grs*mesh->D[nx+mx*mesh->Nq]*mesh->D[my+ny*mesh->Nq];

          id = nx+my*mesh->Nq;
          dfloat Gsr = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G01ID*mesh->Np];
          A[iid] += Gsr*mesh->D[mx+nx*mesh->Nq]*mesh->D[ny+my*mesh->Nq];

          if (nx==mx) {
            for (int k=0;k<mesh->Nq;k++) {
              id = nx+k*mesh->Nq;
              dfloat Gss = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + G11ID*mesh->Np];

              A[iid] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[my+k*mesh->Nq];
            }
          }

          if ((nx==mx)&&(ny==my)) {
            id = nx + ny*mesh->Nq;
            dfloat JW = mesh->ggeo[eM*mesh->Np*mesh->Nggeo + id + GWJID*mesh->Np];
            A[iid] += JW*lambda;
          }
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

