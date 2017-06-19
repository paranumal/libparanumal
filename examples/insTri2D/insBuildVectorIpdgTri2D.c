
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

  dfloat *qmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *qmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *dldxM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *dldxP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *dldyM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *dldyP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint m=0;m<mesh->Nfp;m++) {
	dfloat MSnm = 0;

	for (iint i=0;i<mesh->Np;i++)
	  MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];

	MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]  = MSnm;
      }
    }
  }

  // DrT*MS, DsT*MS
  dfloat *DrTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DsTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (iint f=0;f<mesh->Nfaces;f++) {
    for (iint n=0;n<mesh->Np;n++) {
      for (iint i=0;i<mesh->Nfp;i++) {
        DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        for (iint m=0;m<mesh->Np;m++) {
          DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Dr[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Ds[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
        }
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

    for (iint m=0;m<mesh->Np;m++) {
      for (iint fM=0;fM<mesh->Nfaces;fM++) {
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat hinv = mesh->sgeo[sid+IHID];

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
	
        // extract trace nodes
        for (iint i=0;i<mesh->Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          iint idM    = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
          iint vidM   = mesh->faceNodes[i+fM*mesh->Nfp];
          iint vidP   = mesh->vmapP[idM]%mesh->Np; // only use this to identify location of positive trace vgeo

	  // this needs work
          qmM[i] =0;
          if (vidM == m) qmM[i] =1;
          qmP[i] =0;
          if (vidP == m) qmP[i] =1;
	  
          dldxM[i] = drdx*mesh->Dr[m+vidM*mesh->Np]+dsdx*mesh->Ds[m+vidM*mesh->Np];
	  dldyM[i] = drdy*mesh->Dr[m+vidM*mesh->Np]+dsdy*mesh->Ds[m+vidM*mesh->Np];
          dldxP[i] = drdxP*mesh->Dr[m+vidP*mesh->Np]+dsdxP*mesh->Ds[m+vidP*mesh->Np];
	  dldyP[i] = drdyP*mesh->Dr[m+vidP*mesh->Np]+dsdyP*mesh->Ds[m+vidP*mesh->Np];
        }

        dfloat penalty = tau*hinv; 
        eP = mesh->EToE[eM*mesh->Nfaces+fM];
	
        for (iint n=0;n<mesh->Np;n++) {
	  dfloat SPnm = 0, SxxPnm = 0, SxyPnm = 0, SyxPnm = 0, SyyPnm = 0;
          for (iint i=0;i<mesh->Nfp;i++) {
	    dfloat MSni = MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np];

	    dfloat dlndrMSni = DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np];
	    dfloat dlndsMSni = DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np];
	    
	    dfloat dlndxMSni = drdx*dlndrMSni + dsdx*dlndsMSni;
	    dfloat dlndyMSni = drdy*dlndrMSni + dsdy*dlndsMSni;

	    dfloat dlndxPSni = drdxP*dlndrMSni + dsdxP*dlndsMSni;
	    dfloat dlndyPSni = drdyP*dlndrMSni + dsdyP*dlndsMSni;

	    // penalties on jumps of all derivatives
	    SxxM[m+n*mesh->Np] += -0.5*sJ*MSni*nx*((1+gradqSgn)*dldxM[i]-(1-qSgn)*penalty*nx*qmM[i]);
	    SxyM[m+n*mesh->Np] += -0.5*sJ*MSni*nx*((1+gradqSgn)*dldyM[i]-(1-qSgn)*penalty*ny*qmM[i]);
	    SyxM[m+n*mesh->Np] += -0.5*sJ*MSni*ny*((1+gradqSgn)*dldxM[i]-(1-qSgn)*penalty*nx*qmM[i]);
	    SyyM[m+n*mesh->Np] += -0.5*sJ*MSni*ny*((1+gradqSgn)*dldyM[i]-(1-qSgn)*penalty*ny*qmM[i]);
	    
            SxxM[m+n*mesh->Np] +=  -0.5*(1-qSgn)*sJ*nx*dlndxMSni*qmM[i];
	    SxyM[m+n*mesh->Np] +=  -0.5*(1-qSgn)*sJ*nx*dlndyMSni*qmM[i];
	    SyxM[m+n*mesh->Np] +=  -0.5*(1-qSgn)*sJ*ny*dlndxMSni*qmM[i];
            SyyM[m+n*mesh->Np] +=  -0.5*(1-qSgn)*sJ*ny*dlndyMSni*qmM[i];
	    
	    if(eP>=0) { // when there is a neighbor
	      SxxPnm += -0.5*sJ*MSni*nx*(dldxP[i]+penalty*nx*qmP[i]);
	      SxyPnm += -0.5*sJ*MSni*nx*(dldyP[i]+penalty*ny*qmP[i]);
	      SyxPnm += -0.5*sJ*MSni*ny*(dldxP[i]+penalty*nx*qmP[i]);
	      SyyPnm += -0.5*sJ*MSni*ny*(dldyP[i]+penalty*ny*qmP[i]);
	      
	      SxxPnm +=  +0.5*sJ*nx*dlndxMSni*qmP[i];
	      SyxPnm +=  +0.5*sJ*nx*dlndyMSni*qmP[i];
	      SxyPnm +=  +0.5*sJ*ny*dlndxMSni*qmP[i];
	      SyyPnm +=  +0.5*sJ*ny*dlndyMSni*qmP[i];
	    }
	  }

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

  free(BM);  free(MS);
  free(DrTMS); free(DsTMS);

  free(qmM); free(qmP);
  free(dldxM); free(dldyM);
  free(dldxP); free(dldyP);
}
