
#include "ellipticTri2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

} nonZero_t;

iint addNonZero(nonZero_t *nonZeros, iint nnz, iint row, iint col, iint owner, dfloat val){
  
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
			     iint *BCType, nonZero_t **A, iint *nnzA,
			     hgs_t **hgs, iint *globalStarts, const char *options){

  iint size, rankM;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankM);

  iint Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  iint *globalIds = (iint *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(iint));
  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));

  if (strstr(options,"PROJECT")) {
    // Create a contiguous numbering system, starting from the element-vertex connectivity
    for (iint n=0;n<Nnum;n++) {
      iint id = mesh->gatherLocalIds[n];
      globalIds[id] = mesh->gatherBaseIds[n];
    }

    // squeeze node numbering
    meshParallelConsecutiveGlobalNumbering(Nnum, globalIds, globalOwners, globalStarts);

    //use the ordering to define a gather+scatter for assembly
    *hgs = meshParallelGatherSetup(mesh, Nnum, globalIds, globalOwners);

  } else {
    // every degree of freedom has its own global id
    /* so find number of elements on each rank */
    iint *rankNelements = (iint*) calloc(size, sizeof(iint));
    iint *rankStarts = (iint*) calloc(size+1, sizeof(iint));
    MPI_Allgather(&(mesh->Nelements), 1, MPI_IINT,
		  rankNelements, 1, MPI_IINT, MPI_COMM_WORLD);
    //find offsets
    for(iint r=0;r<size;++r){
      rankStarts[r+1] = rankStarts[r]+rankNelements[r];
    }
    //use the offsets to set a global id
    for (iint e =0;e<mesh->Nelements;e++) {
      for (int n=0;n<mesh->Np;n++) {
        globalIds[e*mesh->Np +n] = n + (e + rankStarts[rankM])*mesh->Np;
        globalOwners[e*mesh->Np +n] = rankM;
      }
    }

    /* do a halo exchange of global node numbers */
    if (mesh->totalHaloPairs) {
      iint *idSendBuffer = (iint *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(iint));
      meshHaloExchange(mesh, mesh->Np*sizeof(iint), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
      free(idSendBuffer);
    }
  }

  iint nnzLocalBound = 4*mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  // drop tolerance for entries in sparse storage

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Np,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Nfp;n++) {
      iint fn = mesh->faceNodes[f*mesh->Nfp+n];

      for (iint m=0;m<mesh->Nfp;m++) {
	iint fm = mesh->faceNodes[f*mesh->Nfp+m];
	dfloat MSnm = 0;

	for (iint i=0;i<mesh->Np;i++){
	  MSnm += mesh->MM[fn+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];
	}
	
	MS[fm+fn*mesh->Np + f*mesh->Np*mesh->Np]  = MSnm;
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
  for(iint eM=0;eM<mesh->Nelements;++eM){

    iint vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(iint n=0;n<mesh->Np;++n){
      for(iint m=0;m<mesh->Np;++m){
        BM[m+n*mesh->Np]  = J*lambda*mesh->MM[m+n*mesh->Np];
        SxxM[m+n*mesh->Np]  = J*drdx*drdx*mesh->Srr[m+n*mesh->Np];
        SxxM[m+n*mesh->Np] += J*drdx*dsdx*mesh->Srs[m+n*mesh->Np];
        SxxM[m+n*mesh->Np] += J*dsdx*drdx*mesh->Ssr[m+n*mesh->Np];
        SxxM[m+n*mesh->Np] += J*dsdx*dsdx*mesh->Sss[m+n*mesh->Np];
			       	      
	SxyM[m+n*mesh->Np]  = J*drdx*drdy*mesh->Srr[m+n*mesh->Np];
        SxyM[m+n*mesh->Np] += J*drdx*dsdy*mesh->Srs[m+n*mesh->Np];
        SxyM[m+n*mesh->Np] += J*dsdx*drdy*mesh->Ssr[m+n*mesh->Np];
        SxyM[m+n*mesh->Np] += J*dsdx*dsdy*mesh->Sss[m+n*mesh->Np];
			       	      
	SyxM[m+n*mesh->Np]  = J*drdy*drdx*mesh->Srr[m+n*mesh->Np];
	SyxM[m+n*mesh->Np] += J*drdy*dsdx*mesh->Srs[m+n*mesh->Np];
        SyxM[m+n*mesh->Np] += J*dsdy*drdx*mesh->Ssr[m+n*mesh->Np];
        SyxM[m+n*mesh->Np] += J*dsdy*dsdx*mesh->Sss[m+n*mesh->Np];
	
	SyyM[m+n*mesh->Np]  = J*drdy*drdy*mesh->Srr[m+n*mesh->Np];
	SyyM[m+n*mesh->Np] += J*drdy*dsdy*mesh->Srs[m+n*mesh->Np];
        SyyM[m+n*mesh->Np] += J*dsdy*drdy*mesh->Ssr[m+n*mesh->Np];
        SyyM[m+n*mesh->Np] += J*dsdy*dsdy*mesh->Sss[m+n*mesh->Np];
      }
    }

    for (iint fM=0;fM<mesh->Nfaces;fM++) {
      // load surface geofactors for this face
      iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];
      dfloat penalty = tau*hinv; 
      
      iint eP = mesh->EToE[eM*mesh->Nfaces+fM];
      if (eP < 0) eP = eM;
      
      iint vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
      
      int qSgn, gradqSgn;
      int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
      iint bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==0){
	qSgn     = 1;
	gradqSgn = 1;
      }else if(bcType==1){ // Dirichlet
	qSgn     = -1;
	gradqSgn =  1;
      } else if (bcType==2){ // Neumann
	qSgn     =  1;
	gradqSgn = -1;
      } else { // Neumann for now
	qSgn     =  1;
	gradqSgn = -1;
      }

      // reset eP
      eP = mesh->EToE[eM*mesh->Nfaces+fM];

      // mass matrix for this face
      dfloat *MSf = MS+fM*mesh->Np*mesh->Np;

      for(iint n=0;n<mesh->Np;++n){
	for(iint m=0;m<mesh->Np;++m){

	  // accumulate face contributions
	  dfloat SPnm = 0, SxxPnm = 0, SxyPnm = 0, SyxPnm = 0, SyyPnm = 0;

	  // OP11 = OP11 + 0.5*( gtau*mmE )
	  dfloat MSfnm = sJ*MSf[n*mesh->Np+m];
	  SxxM[m+n*mesh->Np] += 0.5*MSfnm*nx*nx*(1-qSgn)*penalty;
	  SxyM[m+n*mesh->Np] += 0.5*MSfnm*nx*ny*(1-qSgn)*penalty;
	  SyxM[m+n*mesh->Np] += 0.5*MSfnm*ny*nx*(1-qSgn)*penalty;
	  SyyM[m+n*mesh->Np] += 0.5*MSfnm*ny*ny*(1-qSgn)*penalty;

	  if(eP>=0){
	    // OP12(:,Fm2) = - 0.5*( gtau*mmE(:,Fm1) );
	    SxxPnm += -0.5*MSfnm*nx*nx*penalty;
	    SxyPnm += -0.5*MSfnm*nx*ny*penalty;
	    SyxPnm += -0.5*MSfnm*ny*nx*penalty;
	    SyyPnm += -0.5*MSfnm*ny*ny*penalty;
	  }

	  // now add differential surface terms
	  for(iint i=0;i<mesh->Nfp;++i){
	    iint iM = mesh->faceNodes[fM*mesh->Nfp+i];
	    
	    dfloat MSfni = sJ*MSf[n*mesh->Np+iM]; // surface Jacobian built in
	    dfloat MSfim = sJ*MSf[iM*mesh->Np+m];

	    dfloat DxMim = drdx*mesh->Dr[iM*mesh->Np+m] + dsdx*mesh->Ds[iM*mesh->Np+m];
	    dfloat DyMim = drdy*mesh->Dr[iM*mesh->Np+m] + dsdy*mesh->Ds[iM*mesh->Np+m];
	    dfloat DxMin = drdx*mesh->Dr[iM*mesh->Np+n] + dsdx*mesh->Ds[iM*mesh->Np+n];
	    dfloat DyMin = drdy*mesh->Dr[iM*mesh->Np+n] + dsdy*mesh->Ds[iM*mesh->Np+n];

	    // OP11 = OP11 + 0.5*( - mmE*Dn1)	    
	    SxxM[m+n*mesh->Np] += -0.5*MSfni*nx*(1+gradqSgn)*DxMim;
	    SxyM[m+n*mesh->Np] += -0.5*MSfni*nx*(1+gradqSgn)*DyMim;
	    SyxM[m+n*mesh->Np] += -0.5*MSfni*ny*(1+gradqSgn)*DxMim;
	    SyyM[m+n*mesh->Np] += -0.5*MSfni*ny*(1+gradqSgn)*DyMim;

	    // OP11 = OP11 + (- Dn1'*mmE );
	    SxxM[m+n*mesh->Np] +=  -0.5*nx*(1-qSgn)*DxMin*MSfim;
	    SxyM[m+n*mesh->Np] +=  -0.5*nx*(1-qSgn)*DyMin*MSfim;
	    SyxM[m+n*mesh->Np] +=  -0.5*ny*(1-qSgn)*DxMin*MSfim;
	    SyyM[m+n*mesh->Np] +=  -0.5*ny*(1-qSgn)*DyMin*MSfim;	    
	    
	    if(eP>=0){
	      iint idM = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+n;
	      iint iP  = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo
	      
	      dfloat DxPim = drdxP*mesh->Dr[iP*mesh->Np+m] + dsdxP*mesh->Ds[iP*mesh->Np+m];
	      dfloat DyPim = drdyP*mesh->Dr[iP*mesh->Np+m] + dsdyP*mesh->Ds[iP*mesh->Np+m];
	      
	      //OP12(Fm1,:) = OP12(Fm1,:) - 0.5*(      mmE(Fm1,Fm1)*Dn2(Fm2,:) );
	      SxxPnm += -0.5*MSfni*nx*DxPim;
	      SxyPnm += -0.5*MSfni*nx*DyPim;
	      SyxPnm += -0.5*MSfni*ny*DxPim;
	      SyyPnm += -0.5*MSfni*ny*DyPim;
	      
	      //OP12(:,Fm2) = OP12(:,Fm2) - 0.5*(-Dn1'*mmE(:, Fm1) );
	      SxxPnm +=  +0.5*nx*DxMin*MSfim;
	      SyxPnm +=  +0.5*nx*DyMin*MSfim;
	      SxyPnm +=  +0.5*ny*DxMin*MSfim;
	      SyyPnm +=  +0.5*ny*DxMin*MSfim;
	    }
	  }

	  // store non-zeros for off diagonal block
	  if(eP>=0){
	    // completed  positive trace for this node
	    dfloat S11nmP = SxxPnm + SyyPnm + sigma*(SxxPnm);
	    dfloat S12nmP = sigma*(SxyPnm);
	    dfloat S21nmP = sigma*(SyxPnm);
	    dfloat S22nmP = SxxPnm + SyyPnm + sigma*(SyyPnm);
	    
	    iint row = 2*globalIds[n + eP*mesh->Np];
	    iint col = 2*globalIds[m + eP*mesh->Np];
	    iint owner = globalOwners[n + eP*mesh->Np];
	    
	    nnz = addNonZero(sendNonZeros, nnz, row+0, col+0, owner, S11nmP);
	    nnz = addNonZero(sendNonZeros, nnz, row+1, col+0, owner, S21nmP);
	    nnz = addNonZero(sendNonZeros, nnz, row+0, col+1, owner, S12nmP);
	    nnz = addNonZero(sendNonZeros, nnz, row+1, col+1, owner, S22nmP);
	    
	  }
	}
      }
    }
    
    // store non-zeros for diagonal block
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Np;m++) {
	iint id = m+n*mesh->Np;
	dfloat S11nm = BM[id] + SxxM[id] + SyyM[id] + sigma*SxxM[id]; // diagonal block + laplacian + penalty
	dfloat S12nm = sigma*(SxyM[id]);
	dfloat S21nm = sigma*(SyxM[id]);
	dfloat S22nm = BM[id] + SxxM[id] + SyyM[id] + sigma*SyyM[id]; // diagonal block + laplacian + penalty

	iint row = 2*globalIds[n + eM*mesh->Np];
	iint col = 2*globalIds[m + eM*mesh->Np];
	iint owner = globalOwners[n + eM*mesh->Np];

	nnz = addNonZero(sendNonZeros, nnz, row+0, col+0, owner, S11nm);
	nnz = addNonZero(sendNonZeros, nnz, row+1, col+0, owner, S21nm);
	nnz = addNonZero(sendNonZeros, nnz, row+0, col+1, owner, S12nm);
	nnz = addNonZero(sendNonZeros, nnz, row+1, col+1, owner, S22nm);
      }
    }
  }

  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));
  
  // count how many non-zeros to send to each process
  for(iint n=0;n<nnz;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, nnz, sizeof(nonZero_t), parallelCompareRowColumn);

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_IINT, ArecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  *nnzA = 0;
  for(iint r=0;r<size;++r){
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
  for(iint n=1;n<*nnzA;++n){
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
