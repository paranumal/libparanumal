#include "ellipticQuad2D.h"

// compare on global indices
int parallelCompareRowColumn(const void *a, const void *b){

  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;

  return 0;

}

void ellipticBuildCoarseIpdgQuad2D(mesh2D *mesh, dfloat tau, dfloat lambda, int *BCType, nonZero_t **A, 
                                int *nnzA, int *globalStarts, const char *options);

void ellipticBuildCoarseContinuousQuad2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, int *nnz,
                              hgs_t **hgs, int *globalStarts, const char* options);


void ellipticCoarsePreconditionerSetupQuad2D(mesh_t *mesh, precon_t *precon, dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **V1, nonZero_t **A, int *nnzA,
                                   hgs_t **hgs, int *globalStarts, const char *options){

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int Nnum = mesh->Nverts*(mesh->Nelements);
 
 
  // 2. Build coarse grid element basis functions
  dfloat *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  for(int j=0;j<mesh->Nq;++j){
    for(int i=0;i<mesh->Nq;++i){
      int n = i+j*mesh->Nq;

      dfloat rn = mesh->gllz[i];
      dfloat sn = mesh->gllz[j];

      V1[0*mesh->Np+n] = 0.25*(1-rn)*(1-sn);
      V1[1*mesh->Np+n] = 0.25*(1+rn)*(1-sn);
      V1[2*mesh->Np+n] = 0.25*(1+rn)*(1+sn);
      V1[3*mesh->Np+n] = 0.25*(1-rn)*(1+sn);
    }
  }

  //build coarse grid A
  if (strstr(options,"IPDG")) {
    MPI_Allgather(&(mesh->Nelements), 1, MPI_INT, globalStarts+1, 1, MPI_INT, MPI_COMM_WORLD);
    for(int r=0;r<size;++r)
      globalStarts[r+1] = globalStarts[r]+globalStarts[r+1]*mesh->Nverts;

    ellipticBuildCoarseIpdgQuad2D(mesh,tau,lambda,BCType,&A,&nnz,globalStarts,options);

  } else if (strstr(options,"CONTINUOUS")) {
    ellipticBuildCoarseContinuousQuad2D(mesh,lambda,&A,&nnz,&hgs,globalStarts,options);
  }
}

void ellipticBuildCoarseContinuousQuad2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, int *nnz,
                              hgs_t **hgs, int *globalStarts, const char* options) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ------------------------------------------------------------------------------------
  // 1. Create a contiguous numbering system, starting from the element-vertex connectivity
  int Nnum = mesh->Nverts*mesh->Nelements;
  int *globalNumbering = (int*) calloc(Nnum, sizeof(int));
  int *globalOwners = (int*) calloc(Nnum, sizeof(int));

  // use original vertex numbering
  memcpy(globalNumbering, mesh->EToV, Nnum*sizeof(int));

  // squeeze numbering
  meshParallelConsecutiveGlobalNumbering(Nnum, globalNumbering, globalOwners, globalStarts);

  //use the ordering to define a gather+scatter for assembly
  *hgs = meshParallelGatherSetup(mesh, Nnum, globalNumbering, globalOwners);

  // ------------------------------------------------------------------------------------
  // Build non-zeros of stiffness matrix (unassembled)
  int nnzLocal = mesh->Nverts*mesh->Nverts*(mesh->Nelements+mesh->totalHaloPairs);
  int   *rowsA = (int*) calloc(nnzLocal, sizeof(int));
  int   *colsA = (int*) calloc(nnzLocal, sizeof(int));
  dfloat *valsA = (dfloat*) calloc(nnzLocal, sizeof(dfloat));

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocal, sizeof(nonZero_t));
  int *AsendCounts  = (int*) calloc(size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(size, sizeof(int));
  int *AsendOffsets = (int*) calloc(size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(size+1, sizeof(int));

  dfloat *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vr1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vs1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  for(int j=0;j<mesh->Nq;++j){
    for(int i=0;i<mesh->Nq;++i){
      int n = i+j*mesh->Nq;

      dfloat rn = mesh->gllz[i];
      dfloat sn = mesh->gllz[j];
      V1[0*mesh->Np+n] = 0.25*(1-rn)*(1-sn);
      V1[1*mesh->Np+n] = 0.25*(1+rn)*(1-sn);
      V1[2*mesh->Np+n] = 0.25*(1+rn)*(1+sn);
      V1[3*mesh->Np+n] = 0.25*(1-rn)*(1+sn);

      Vr1[0*mesh->Np+n] = 0.25*(-1)*(1-sn);
      Vr1[1*mesh->Np+n] = 0.25*(+1)*(1-sn);
      Vr1[2*mesh->Np+n] = 0.25*(+1)*(1+sn);
      Vr1[3*mesh->Np+n] = 0.25*(-1)*(1+sn);

      Vs1[0*mesh->Np+n] = 0.25*(1-rn)*(-1);
      Vs1[1*mesh->Np+n] = 0.25*(1+rn)*(-1);
      Vs1[2*mesh->Np+n] = 0.25*(1+rn)*(+1);
      Vs1[3*mesh->Np+n] = 0.25*(1-rn)*(+1);
    }
  }

  int cnt = 0;
  printf("Building coarse matrix system\n");
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Nverts;++n){
      for(int m=0;m<mesh->Nverts;++m){
      	dfloat Snm = 0;

      	// use GLL nodes for integration
      	// (since Jacobian is high order tensor-product polynomial)
      	for(int j=0;j<mesh->Nq;++j){
      	  for(int i=0;i<mesh->Nq;++i){
      	    int id = i+j*mesh->Nq;

      	    dfloat Vr1ni = Vr1[n*mesh->Np+id];
      	    dfloat Vs1ni = Vs1[n*mesh->Np+id];
      	    dfloat V1ni  = V1[n*mesh->Np+id];

      	    dfloat Vr1mi = Vr1[m*mesh->Np+id];
      	    dfloat Vs1mi = Vs1[m*mesh->Np+id];
      	    dfloat V1mi  = V1[m*mesh->Np+id];

      	    dfloat rx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + RXID*mesh->Np];
      	    dfloat sx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + SXID*mesh->Np];
      	    dfloat ry = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + RYID*mesh->Np];
      	    dfloat sy = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + SYID*mesh->Np];
      	    dfloat JW = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + JWID*mesh->Np];

      	    dfloat Vx1ni = rx*Vr1ni+sx*Vs1ni;
      	    dfloat Vy1ni = ry*Vr1ni+sy*Vs1ni;
      	    dfloat Vx1mi = rx*Vr1mi+sx*Vs1mi;
      	    dfloat Vy1mi = ry*Vr1mi+sy*Vs1mi;

      	    Snm += (Vx1ni*Vx1mi+Vy1ni*Vy1mi)*JW;
      	    Snm += (lambda*V1ni*V1mi)*JW;
      	  }
      	}

        dfloat nonZeroThreshold = 1e-7;
        if(fabs(Snm)>nonZeroThreshold) {
          valsA[cnt] = Snm;
          rowsA[cnt] = e*mesh->Nverts+n;
          colsA[cnt] = e*mesh->Nverts+m;

        	// pack non-zero
        	sendNonZeros[cnt].val = Snm;
        	sendNonZeros[cnt].row = globalNumbering[e*mesh->Nverts+n];
        	sendNonZeros[cnt].col = globalNumbering[e*mesh->Nverts+m];
        	sendNonZeros[cnt].ownerRank = globalOwners[e*mesh->Nverts+n];

        	++cnt;
        }
      }
    }
  }

  // count how many non-zeros to send to each process
  for(int n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_INT, ArecvCounts, 1, MPI_INT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnz = 0;
  for(int r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnz += ArecvCounts[r]/sizeof(nonZero_t);
  }

  *A = (nonZero_t*) calloc(*nnz, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
		(*A), ArecvCounts, ArecvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(int n=1;n<*nnz;++n){
    if((*A)[n].row == (*A)[cnt].row &&
       (*A)[n].col == (*A)[cnt].col){
      (*A)[cnt].val += (*A)[n].val;
    }
    else{
      ++cnt;
      (*A)[cnt] = (*A)[n];
    }
  }
  *nnz = cnt+1;
  printf("Done building coarse matrix system\n");

  free(globalNumbering); free(globalOwners);
  free(V1); free(Vr1); free(Vs1);
  free(rowsA); free(colsA); free(valsA);
  free(sendNonZeros);

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);
}

void ellipticBuildCoarseIpdgQuad2D(mesh2D *mesh, dfloat tau, dfloat lambda, int *BCType,
                      nonZero_t **A, int *nnzA, int *globalStarts, const char *options){

  int size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  int Nnum = mesh->Nverts*mesh->Nelements;

  // create a global numbering system
  int *globalIds = (int *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts,sizeof(int));

  // every degree of freedom has its own global id
  /* so find number of elements on each rank */
  int *rankNelements = (int*) calloc(size, sizeof(int));
  int *rankStarts = (int*) calloc(size+1, sizeof(int));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_INT,
    rankNelements, 1, MPI_INT, MPI_COMM_WORLD);
  //find offsets
  for(int r=0;r<size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (int e =0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Nverts;n++) {
      globalIds[e*mesh->Nverts +n] = n + (e + rankStarts[rankM])*mesh->Nverts;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    int *idSendBuffer = (int *) calloc(mesh->Nverts*mesh->totalHaloPairs,sizeof(int));
    meshHaloExchange(mesh, mesh->Nverts*sizeof(int), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Nverts);
    free(idSendBuffer);
  }

  // Build coarse basis restriction
  dfloat *V  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  for(int j=0;j<mesh->Nq;++j){
    for(int i=0;i<mesh->Nq;++i){
      int n = i+j*mesh->Nq;

      dfloat rn = mesh->gllz[i];
      dfloat sn = mesh->gllz[j];
      V[0*mesh->Np+n] = 0.25*(1-rn)*(1-sn);
      V[1*mesh->Np+n] = 0.25*(1+rn)*(1-sn);
      V[2*mesh->Np+n] = 0.25*(1+rn)*(1+sn);
      V[3*mesh->Np+n] = 0.25*(1-rn)*(1+sn);
    }
  }

  int nnzLocalBound = mesh->Nverts*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-12;

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  dfloat *BP = (dfloat *) calloc(mesh->Np*mesh->Np*mesh->Nfaces,sizeof(dfloat));

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
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

  // reset non-zero counter
  int nnz = 0;

  // loop over all elements
  for(int eM=0;eM<mesh->Nelements;++eM){

    //zero out BP
    for (int fM=0;fM<mesh->Nfaces;fM++)
      for (int n=0;n<mesh->Np;n++)
        for (int m=0;m<mesh->Np;m++)
          BP[m+n*mesh->Np+fM*mesh->Np*mesh->Np] = 0;

    /* build Dx,Dy (forget the TP for the moment) */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh->Np;++i){
          int base = eM*mesh->Np*mesh->Nvgeo + i;
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
          BM[m+n*mesh->Np] = JW*(dlndx*dlmdx+dlndy*dlmdy);
          BM[m+n*mesh->Np] += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh->Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh->Nfp;++i){
            int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

            // grab vol geofacs at surface nodes
            int baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
            dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
            dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
            dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
            dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];

            // double check vol geometric factors are in halo storage of vgeo
            int idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
            int vidP    = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo
            int localEP = mesh->vmapP[idM]/mesh->Np;
            int baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
            dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
            dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
            dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];

            // grab surface geometric factors
            int base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
            dfloat nx = mesh->sgeo[base+NXID];
            dfloat ny = mesh->sgeo[base+NYID];
            dfloat wsJ = mesh->sgeo[base+WSJID];
            dfloat hinv = mesh->sgeo[base+IHID];

            // form negative trace terms in IPDG
            int idnM = n*mesh->Np+vidM;
            int idmM = m*mesh->Np+vidM;
            int idmP = m*mesh->Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM;
            dfloat lmM = B[idmM];

            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP];
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP];
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP;
            dfloat lmP = B[idmP];

            dfloat penalty = tau*(mesh->N+1)*(mesh->N+1)*hinv;

            BM[m+n*mesh->Np] += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            BM[m+n*mesh->Np] += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            BM[m+n*mesh->Np] += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            int eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            if (eP < 0) {
              int qSgn, gradqSgn;
              int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
              int bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
              if(bcType==1){ // Dirichlet
                qSgn     = -1;
                gradqSgn =  1;
              } else if (bcType==2){ // Neumann
                qSgn     =  1;
                gradqSgn = -1;
              } else { // Neumann for now
                qSgn     =  1;
                gradqSgn = -1;
              }

              BM[m+n*mesh->Np] += -0.5*gradqSgn*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
              BM[m+n*mesh->Np] += +0.5*qSgn*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
              BM[m+n*mesh->Np] += -0.5*qSgn*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-)
            } else {
              AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
              AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
              AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
            BP[m+n*mesh->Np+fM*mesh->Np*mesh->Np] = AnmP;
          }
        }
      }
    }

    //transform to coarse operator Ac = VAV^T
    for (int fM=0;fM<mesh->Nfaces;fM++) {
      int eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) eP = eM;

      for (int j=0;j<mesh->Nverts;j++) {
        for (int i=0;i<mesh->Nverts;i++) {
          dfloat AnmP = 0;
          for (int m=0;m<mesh->Np;m++) {
            for (int n=0;n<mesh->Np;n++) {
              AnmP += V[n+i*mesh->Np]*BP[m+n*mesh->Np+fM*mesh->Np*mesh->Np]*V[m+j*mesh->Np];
            }
          }
          if(fabs(AnmP)>tol){
            (*A)[nnz].row = globalIds[eM*mesh->Nverts + i];
            (*A)[nnz].col = globalIds[eP*mesh->Nverts + j];
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = globalOwners[eM*mesh->Nverts + i];
            ++nnz;
          }
        }
      }
    }

    for (int j=0;j<mesh->Nverts;j++) {
      for (int i=0;i<mesh->Nverts;i++) {
        dfloat Anm = 0;
        for (int m=0;m<mesh->Np;m++) {
          for (int n=0;n<mesh->Np;n++) {
            Anm += V[n+i*mesh->Np]*BM[m+n*mesh->Np]*V[m+j*mesh->Np];
          }
        }

        if(fabs(Anm)>tol){
          (*A)[nnz].row = globalIds[eM*mesh->Nverts+i];
          (*A)[nnz].col = globalIds[eM*mesh->Nverts+j];
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = globalOwners[eM*mesh->Nverts + i];
          ++nnz;
        }
      }
    }
  }

  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnzA, sizeof(nonZero_t), parallelCompareRowColumn);

  free(globalIds);
  free(V); free(BP); free(BM);
  free(B);  free(Br); free(Bs);
}