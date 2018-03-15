
#include "ins2D.h"

int addNonZero(nonZero_t *nonZeros, int nnz, int row, int col, int owner, dfloat val){
  
  dfloat tol = 1e-12;
  
  if(fabs(val)>tol){
    nonZeros[nnz].val = val;
    nonZeros[nnz].row = row;
    nonZeros[nnz].col = col;
    nonZeros[nnz].ownerRank = owner;
    ++nnz;
  }
  
  return nnz;
}


int parallelCompareRowColumn(const void *a, const void *b);

void insBuildVectorIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat sigma, dfloat lambda,
			     int *BCType, nonZero_t **A, int *nnzA,
			     hgs_t **hgs, int *globalStarts, const char *options){

  int size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  int Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  int *globalIds = (int *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(int));
  int *globalOwners = (int*) calloc(Nnum, sizeof(int));

  if (strstr(options,"PROJECT")) {
    // Create a contiguous numbering system, starting from the element-vertex connectivity
    for (int n=0;n<Nnum;n++) {
      int id = mesh->gatherLocalIds[n];
      globalIds[id] = mesh->gatherBaseIds[n];
    }

    // squeeze node numbering
    meshParallelConsecutiveGlobalNumbering(Nnum, globalIds, globalOwners, globalStarts);

    //use the ordering to define a gather+scatter for assembly
    *hgs = meshParallelGatherSetup(mesh, Nnum, globalIds, globalOwners);

  } else {
    // every degree of freedom has its own global id
    /* so find number of elements on each rank */
    int *rankNelements = (int*) calloc(size, sizeof(int));
    int *rankStarts = (int*) calloc(size+1, sizeof(int));
    MPI_Allgather(&(mesh->Nelements), 1, MPI_int,
		  rankNelements, 1, MPI_int, MPI_COMM_WORLD);
    //find offsets
    for(int r=0;r<size;++r){
      rankStarts[r+1] = rankStarts[r]+rankNelements[r];
    }
    //use the offsets to set a global id
    for (int e =0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Np;n++) {
        globalIds[e*mesh->Np +n] = n + (e + rankStarts[rankM])*mesh->Np;
        globalOwners[e*mesh->Np +n] = rankM;
      }
    }

    /* do a halo exchange of global node numbers */
    if (mesh->totalHaloPairs) {
      int *idSendBuffer = (int *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(int));
      meshHaloExchange(mesh, mesh->Np*sizeof(int), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
      free(idSendBuffer);
    }
  }

  int nnzLocalBound = 4*mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  // drop tolerance for entries in sparse storage

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Nfp*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Nfp;n++) {
      int fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (int m=0;m<mesh->Nfp;m++) {
	int fm = mesh->faceNodes[f*mesh->Nfp+m];
	dfloat MSnm = 0;

	for (int i=0;i<mesh->Np;i++){
	  MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
	}
	
	MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Nfp]  = MSnm;
      }
    }
  }

  
  // reset non-zero counter
  int nnz = 0;

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  dfloat *SxxM = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *SxyM = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *SyxM = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *SyyM = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  // loop over all elements
  for(int eM=0;eM<mesh->Nelements;++eM){

    int vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        BM[n*mesh->Np+m]  = J*lambda*mesh->MM[n*mesh->Np+m];
        SxxM[n*mesh->Np+m]  = J*drdx*drdx*mesh->Srr[n*mesh->Np+m];
        SxxM[n*mesh->Np+m] += J*drdx*dsdx*mesh->Srs[n*mesh->Np+m];
        SxxM[n*mesh->Np+m] += J*dsdx*drdx*mesh->Ssr[n*mesh->Np+m];
        SxxM[n*mesh->Np+m] += J*dsdx*dsdx*mesh->Sss[n*mesh->Np+m];
			       	      
	SxyM[n*mesh->Np+m]  = J*drdx*drdy*mesh->Srr[n*mesh->Np+m];
        SxyM[n*mesh->Np+m] += J*drdx*dsdy*mesh->Srs[n*mesh->Np+m];
        SxyM[n*mesh->Np+m] += J*dsdx*drdy*mesh->Ssr[n*mesh->Np+m];
        SxyM[n*mesh->Np+m] += J*dsdx*dsdy*mesh->Sss[n*mesh->Np+m];
			       	      
	SyxM[n*mesh->Np+m]  = J*drdy*drdx*mesh->Srr[n*mesh->Np+m];
	SyxM[n*mesh->Np+m] += J*drdy*dsdx*mesh->Srs[n*mesh->Np+m];
        SyxM[n*mesh->Np+m] += J*dsdy*drdx*mesh->Ssr[n*mesh->Np+m];
        SyxM[n*mesh->Np+m] += J*dsdy*dsdx*mesh->Sss[n*mesh->Np+m];
	
	SyyM[n*mesh->Np+m]  = J*drdy*drdy*mesh->Srr[n*mesh->Np+m];
	SyyM[n*mesh->Np+m] += J*drdy*dsdy*mesh->Srs[n*mesh->Np+m];
        SyyM[n*mesh->Np+m] += J*dsdy*drdy*mesh->Ssr[n*mesh->Np+m];
        SyyM[n*mesh->Np+m] += J*dsdy*dsdy*mesh->Sss[n*mesh->Np+m];
      }
    }

   
    for (int fM=0;fM<mesh->Nfaces;fM++) {
      
      dfloat *SxxP = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
      dfloat *SxyP = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
      dfloat *SyxP = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
      dfloat *SyyP = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
      
      // load surface geofactors for this face
      int sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];
      dfloat penalty = tau*hinv; 
      
      int eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) eP = eM;
      
      int vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
      
      int bcD, bcN;
      int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
      int bcType = 0;

      if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType<=0){
	bcD = 0;
	bcN = 0;
      }else if(bcType==1){ // Dirichlet
	bcD = 1;
	bcN = -1;
      } else if (bcType==2){ // Neumann
	bcD = -1;
	bcN = 1;
      } else { // Neumann for now
	bcD = -1;
	bcN = 1;
      }
      
      // reset eP
      eP = mesh->EToE[eM*mesh->Nfaces+fM];

      // mass matrix for this face
      dfloat *MSf = MS+fM*mesh->Nfp*mesh->Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<mesh->Nfp;++n){
	for(int m=0;m<mesh->Nfp;++m){
	  int nM = mesh->faceNodes[fM*mesh->Nfp+n];
	  int mM = mesh->faceNodes[fM*mesh->Nfp+m];

	  // OP11 = OP11 + 0.5*( gtau*mmE )
	  dfloat MSfnm = sJ*MSf[n*mesh->Nfp+m];
	  SxxM[nM*mesh->Np+mM] += 0.5*nx*nx*(1+bcD)*penalty*MSfnm;
	  SxyM[nM*mesh->Np+mM] += 0.5*nx*ny*(1+bcD)*penalty*MSfnm;
	  SyxM[nM*mesh->Np+mM] += 0.5*ny*nx*(1+bcD)*penalty*MSfnm;
	  SyyM[nM*mesh->Np+mM] += 0.5*ny*ny*(1+bcD)*penalty*MSfnm;

	  // neighbor penalty term
	  if(eP>=0){
	    int idM = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+n;
	    int mP  = mesh->vmapP[idM]%mesh->Np; 

	    // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
	    SxxP[nM*mesh->Np+mP] += -0.5*nx*nx*penalty*MSfnm;
	    SxyP[nM*mesh->Np+mP] += -0.5*nx*ny*penalty*MSfnm;
	    SyxP[nM*mesh->Np+mP] += -0.5*ny*nx*penalty*MSfnm;
	    SyyP[nM*mesh->Np+mP] += -0.5*ny*ny*penalty*MSfnm;
	  }
	}
      }
    
      // now add differential surface terms
      for(int n=0;n<mesh->Nfp;++n){
	for(int m=0;m<mesh->Np;++m){
	  int nM = mesh->faceNodes[fM*mesh->Nfp+n];
	  
	  for(int i=0;i<mesh->Nfp;++i){
	    int iM = mesh->faceNodes[fM*mesh->Nfp+i];
	    int iP = mesh->vmapP[i + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;
	      
	    dfloat MSfni = sJ*MSf[n*mesh->Nfp+i]; // surface Jacobian built in
	    
	    dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+m] + dsdx*mesh->Ds[iM*mesh->Np+m];
	    dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+m] + dsdy*mesh->Ds[iM*mesh->Np+m];
	    dfloat DxPim = drdxP*mesh->Dr[iP*mesh->Np+m] + dsdxP*mesh->Ds[iP*mesh->Np+m];
	    dfloat DyPim = drdyP*mesh->Dr[iP*mesh->Np+m] + dsdyP*mesh->Ds[iP*mesh->Np+m];

	    // OP11 = OP11 + 0.5*( - mmE*Dn1)	    
	    SxxM[nM*mesh->Np+m] += -0.5*nx*(1+bcN)*MSfni*DxMim;
	    SxyM[nM*mesh->Np+m] += -0.5*nx*(1+bcN)*MSfni*DyMim;
	    SyxM[nM*mesh->Np+m] += -0.5*ny*(1+bcN)*MSfni*DxMim;
	    SyyM[nM*mesh->Np+m] += -0.5*ny*(1+bcN)*MSfni*DyMim;
	    
	    if(eP>=0){
	      //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
	      SxxP[nM*mesh->Np+m] += -0.5*nx*MSfni*DxPim;
	      SxyP[nM*mesh->Np+m] += -0.5*nx*MSfni*DyPim;
	      SyxP[nM*mesh->Np+m] += -0.5*ny*MSfni*DxPim;
	      SyyP[nM*mesh->Np+m] += -0.5*ny*MSfni*DyPim;
	    }
	  }
	}
      }
    
      for(int n=0;n<mesh->Np;++n){
	for(int m=0;m<mesh->Nfp;++m){
	  int mM = mesh->faceNodes[fM*mesh->Nfp+m];
	  int mP = mesh->vmapP[m + fM*mesh->Nfp+eM*mesh->Nfp*mesh->Nfaces]%mesh->Np;
	  
	  for(int i=0;i<mesh->Nfp;++i){
	    int iM = mesh->faceNodes[fM*mesh->Nfp+i];	

	    dfloat MSfim = sJ*MSf[i*mesh->Nfp+m];
	    
	    dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
	    dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];
	  
	    // OP11 = OP11 + (- Dn1'*mmE );
	    SxxM[n*mesh->Np+mM] +=  -0.5*nx*(1+bcD)*DxMin*MSfim;
	    SxyM[n*mesh->Np+mM] +=  -0.5*nx*(1+bcD)*DyMin*MSfim;
	    SyxM[n*mesh->Np+mM] +=  -0.5*ny*(1+bcD)*DxMin*MSfim;
	    SyyM[n*mesh->Np+mM] +=  -0.5*ny*(1+bcD)*DyMin*MSfim;

	    if(eP>=0){
	      //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
	      SxxP[n*mesh->Np+mP] +=  +0.5*nx*DxMin*MSfim;
	      SyxP[n*mesh->Np+mP] +=  +0.5*nx*DyMin*MSfim;
	      SxyP[n*mesh->Np+mP] +=  +0.5*ny*DxMin*MSfim;
	      SyyP[n*mesh->Np+mP] +=  +0.5*ny*DyMin*MSfim;
	    }
	  }
	}
      }

      // store non-zeros for off diagonal block
      if(eP>=0){
	for(int n=0;n<mesh->Np;++n){
	  for(int m=0;m<mesh->Np;++m){
	    int id = n*mesh->Np+m;
	    
	    // completed  positive trace for this node
	    dfloat S11nmP = SxxP[id] + SyyP[id] + sigma*(SxxP[id]);
	    dfloat S12nmP = sigma*(SxyP[id]);
	    dfloat S21nmP = sigma*(SyxP[id]);
	    dfloat S22nmP = SxxP[id] + SyyP[id] + sigma*(SyyP[id]);
	    
	    int row = 2*globalIds[n + eP*mesh->Np];
	    int col = 2*globalIds[m + eP*mesh->Np];
	    int owner = globalOwners[n + eP*mesh->Np];
	    
	    nnz = addNonZero(sendNonZeros, nnz, row+0, col+0, owner, S11nmP);
	    nnz = addNonZero(sendNonZeros, nnz, row+1, col+0, owner, S21nmP);
	    nnz = addNonZero(sendNonZeros, nnz, row+0, col+1, owner, S12nmP);
	    nnz = addNonZero(sendNonZeros, nnz, row+1, col+1, owner, S22nmP);
	    
	  }
	}
      }
      
      free(SxxP); free(SxyP); free(SyxP); free(SyyP);
    }
    
    // store non-zeros for diagonal block
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
	int id = n*mesh->Np+m;
	dfloat S11nm = BM[id] + SxxM[id] + SyyM[id] + sigma*SxxM[id]; // diagonal block + laplacian + penalty
	dfloat S12nm = sigma*(SxyM[id]);
	dfloat S21nm = sigma*(SyxM[id]);
	dfloat S22nm = BM[id] + SxxM[id] + SyyM[id] + sigma*SyyM[id]; // diagonal block + laplacian + penalty
	
	int row = 2*globalIds[n + eM*mesh->Np];
	int col = 2*globalIds[m + eM*mesh->Np];
	int owner = globalOwners[n + eM*mesh->Np];
	
	nnz = addNonZero(sendNonZeros, nnz, row+0, col+0, owner, S11nm);
	nnz = addNonZero(sendNonZeros, nnz, row+1, col+0, owner, S21nm);
	nnz = addNonZero(sendNonZeros, nnz, row+0, col+1, owner, S12nm);
	nnz = addNonZero(sendNonZeros, nnz, row+1, col+1, owner, S22nm);
      }
    }
  }
  
  int *AsendCounts  = (int*) calloc(size, sizeof(int));
  int *ArecvCounts  = (int*) calloc(size, sizeof(int));
  int *AsendOffsets = (int*) calloc(size+1, sizeof(int));
  int *ArecvOffsets = (int*) calloc(size+1, sizeof(int));
  
  // count how many non-zeros to send to each process
  for(int n=0;n<nnz;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_int, ArecvCounts, 1, MPI_int, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnzA = 0;
  for(int r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    *nnzA += ArecvCounts[r]/sizeof(nonZero_t);
  }

  *A = (nonZero_t*) calloc(*nnzA, sizeof(nonZero_t));

  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
		(*A), ArecvCounts, ArecvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), *nnzA, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  nnz = 0;
  for(int n=1;n<*nnzA;++n){
    if((*A)[n].row == (*A)[nnz].row &&
       (*A)[n].col == (*A)[nnz].col){
      (*A)[nnz].val += (*A)[n].val;
    }
    else{
      ++nnz;
      (*A)[nnz] = (*A)[n];
    }
  }
  *nnzA = nnz+1;

  free(globalIds);
  free(globalOwners);
  free(sendNonZeros);
  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);

  free(BM);
  free(MS);

}
