#include "ellipticHex3D.h"

void BuildLocalIpdgDiag(mesh3D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat* Bs, dfloat* Bt, dlong eM, dfloat *A);

void BuildLocalContinuousDiag(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat* Bs, dfloat* Bt, dfloat *A);

void ellipticBuildJacobiHex3D(solver_t* solver, dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options){

  mesh3D *mesh = solver->mesh;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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

  dlong diagNnum = mesh->Np*mesh->Nelements;

  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  if(rank==0) printf("Building diagonal...");fflush(stdout);

  // loop over all elements
  #pragma omp parallel
  {
    #pragma omp for 
    for(dlong eM=0;eM<mesh->Nelements;++eM){
      //build the patch A matrix for this element
      if (strstr(options,"IPDG")) {
        BuildLocalIpdgDiag(mesh, tau, lambda, BCType, B, Br, Bs, Bt, eM, diagA + eM*mesh->Np);
      } else if (strstr(options,"CONTINUOUS")) {
        BuildLocalContinuousDiag(solver, lambda, eM, B, Br, Bs, Bt, diagA + eM*mesh->Np);
      }
    }
  }

  if (strstr(options,"CONTINUOUS")) 
    gsParallelGatherScatter(mesh->hostGsh, diagA, dfloatString, "add"); 
    
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<mesh->Nelements*mesh->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(rank==0) printf("done.\n");

  free(B); free(Br); free(Bs); free(Bt);
  free(diagA);
}



//returns the patch A matrix for element eM
void BuildLocalIpdgDiag(mesh3D *mesh, dfloat tau, dfloat lambda, int* BCType,
                        dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dlong eM, dfloat *diagA) {

  /* start with stiffness matrix  */
  for(int n=0;n<mesh->Np;++n){
    diagA[n] = 0;

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
      dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
      dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
      dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];    
      diagA[n] += JW*(dlndx*dlndx+dlndy*dlndy+dlndz*dlndz);
      diagA[n] += lambda*JW*B[idn]*B[idn];
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
        dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
        dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
        dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
        dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
        dfloat lnM = B[idnM];
        
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

        diagA[n] += -0.5*(1+bcD)*(1-bcN)*wsJ*lnM*ndotgradlnM;  // -(ln^-, N.grad lm^-)
        diagA[n] += -0.5*(1+bcD)*(1-bcN)*wsJ*ndotgradlnM*lnM;  // -(N.grad ln^-, lm^-)
        diagA[n] += +0.5*(1+bcD)*(1-bcN)*wsJ*penalty*lnM*lnM; // +((tau/h)*ln^-,lm^-)
      }
    }
  }
}

void BuildLocalContinuousDiag(solver_t* solver, dfloat lambda,
                                  dlong eM, dfloat *B, dfloat *Br, dfloat* Bs, dfloat* Bt, dfloat *diagA) {

  mesh3D *mesh = solver->mesh;

  for (int nz=0;nz<mesh->Nq;nz++) {
  for (int ny=0;ny<mesh->Nq;ny++) {
  for (int nx=0;nx<mesh->Nq;nx++) {
    int idn = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
    if (solver->mapB[idn+eM*mesh->Np]!=1) {            
      diagA[idn] = 0;

      int id = nx+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
      dlong base = eM*mesh->Np*mesh->Nggeo;

    
      dfloat Grs = mesh->ggeo[base + id + G01ID*mesh->Np];
      diagA[idn] += 2*Grs*mesh->D[nx+nx*mesh->Nq]*mesh->D[ny+ny*mesh->Nq];

      dfloat Grt = mesh->ggeo[base + id + G02ID*mesh->Np];
      diagA[idn] += 2*Grt*mesh->D[nx+nx*mesh->Nq]*mesh->D[nz+nz*mesh->Nq];
    
      dfloat Gst = mesh->ggeo[base + id + G12ID*mesh->Np];
      diagA[idn] += 2*Gst*mesh->D[ny+ny*mesh->Nq]*mesh->D[nz+nz*mesh->Nq];

      for (int k=0;k<mesh->Nq;k++) {
        int iid = k+ny*mesh->Nq+nz*mesh->Nq*mesh->Nq;
        dfloat Grr = mesh->ggeo[base + iid + G00ID*mesh->Np];
        diagA[idn] += Grr*mesh->D[nx+k*mesh->Nq]*mesh->D[nx+k*mesh->Nq];
      }

      for (int k=0;k<mesh->Nq;k++) {
        int iid = nx+k*mesh->Nq+nz*mesh->Nq*mesh->Nq;
        dfloat Gss = mesh->ggeo[base + iid + G11ID*mesh->Np];
        diagA[idn] += Gss*mesh->D[ny+k*mesh->Nq]*mesh->D[ny+k*mesh->Nq];
      }
    
      for (int k=0;k<mesh->Nq;k++) {
        int iid = nx+ny*mesh->Nq+k*mesh->Nq*mesh->Nq;
        dfloat Gtt = mesh->ggeo[base + iid + G22ID*mesh->Np];
        diagA[idn] += Gtt*mesh->D[nz+k*mesh->Nq]*mesh->D[nz+k*mesh->Nq];
      }
      
      dfloat JW = mesh->ggeo[base + id + GWJID*mesh->Np];
      diagA[idn] += JW*lambda;
    } else {
      diagA[idn] = 1; //just put a 1 so diagA is invertable
    }
  }
  }
  }

  //add the rank boost for the allNeumann Poisson problem
  if (solver->allNeumann) {
    for(int n=0;n<mesh->Np;++n){
      if (solver->mapB[n+eM*mesh->Np]!=1) { //dont fill rows for masked nodes
        diagA[n] += solver->allNeumannPenalty*solver->allNeumannScale*solver->allNeumannScale;
      } 
    }
  }
}