#include "ellipticQuad2D.h"

int compareRowCol(void *a, void *b){
  
  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->row < fb->row) return -1;
  if(fa->row > fb->row) return +1;

  if(fa->col < fb->col) return -1;
  if(fa->col > fb->col) return +1;
  
  return 0;  
}

int compareOwner(void *a, void *b){
  
  nonZero_t *fa = (nonZero_t*) a;
  nonZero_t *fb = (nonZero_t*) b;

  if(fa->own < fb->own) return -1;
  if(fa->own > fb->own) return +1;
  
  return 0;  
}

void ellipticCoarsePreconditionerSetupQuad2D(mesh_t *mesh, precon_t *precon, dfloat lambda){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // ------------------------------------------------------------------------------------
  // 0. Create a contiguous numbering system, starting from the element-vertex connectivity
  
  iint Nnum = mesh->Nverts*mesh->Nelements;
  
  iint *globalNumbering = (iint*) calloc(Nnum, sizeof(iint));
  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));

  // use original vertex numbering
  memcpy(globalNumbering, mesh->EToV, Nnum);

  // squeeze numbering
  meshParallelConsecutiveGlobalNumbering(Nnum, globalNumbering, globalOwners);

  // ------------------------------------------------------------------------------------
  // 1. Build non-zeros of stiffness matrix (unassembled)

  // counts
  iint *sendCounts = (iint*) calloc(size, sieof(iint));
  iint *recvCounts = (iint*) calloc(size, sieof(iint));
  
  // build stiffness matrix for vertex modes
  iint NnonZeros = mesh->Nverts*mesh->Nverts*mesh->Nelements;
  nonZero_t *sendNonZeros = (nonZero_t*) calloc(NnonZeros, sizeof(nonZero_t));
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      for(iint m=0;m<mesh->Nverts;++m){
	dfloat Snm = 0;

	// use GLL nodes for integration (since Jacobian is high order tensor-product polynomial)
	for(iint j=0;j<mesh->Nq;++j){
	  for(iint i=0;i<mesh->Nq;++i){
	    iint id = i+j*mesh->Nq;
	    
	    dfloat Vr1ni = mesh->Vr1[n*mesh->Np+id];
	    dfloat Vs1ni = mesh->Vs1[n*mesh->Np+id];
	    
	    dfloat Vr1mi = mesh->Vr1[m*mesh->Np+id];
	    dfloat Vs1mi = mesh->Vs1[m*mesh->Np+id];

	    dfloat V1ni  = mesh->V1[n*mesh->Np+id];
	    dfloat V1mi  = mesh->V1[m*mesh->Np+id];

	    dfloat rx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + RXID*mesh->Np];
	    dfloat sx = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + SXID*mesh->Np];
	    dfloat ry = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + RYID*mesh->Np];
	    dfloat sy = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + SYID*mesh->Np];
	    dfloat JW = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + id + JWID*mesh->Np];
	    
	    dfloat Vx1ni = rx*Vr1ni+sx*Vs1ni;
	    dfloat Vy1ni = ry*Vr1ni+sy*Vs1ni;
	    dfloat Vx1mi = rx*Vr1mi+sx*Vs1mi;
	    dfloat Vy1mi = ry*Vr1mi+sy*Vs1mi;
	    
	    Snm += (Vx1ni*Vx1mi+Vy1ni*Vy1mi + lambda*V1ni*V1mi)*JW;
	  }
	}
	
	iint row = globalNumbering[e*mesh->Nverts+n];
	iint col = globalNumbering[e*mesh->Nverts+m];
	iint own = globalOwners[e*mesh->Nverts+m]; // own by row
	sendNonZeros[cnt].row = row;
	sendNonZeros[cnt].col = col;
	sendNonZeros[cnt].val = Snm;
	sendNonZeros[cnt].own = own;

	// increment the send counter for the owner rank
	++sendCounts[own];
	
	++cnt;
      }
    }
  }

  // ------------------------------------------------------------------------------------
  // 2. Sort non-zeros based on owner rank and send to owner
  
  // sort the non-zeros by owner rank
  qsort(sendNonZeros, NnonZeros, sizeof(nonZero_t), compareOwner);
  
  iint *sendOffsets = (iint*) calloc(size+1, sieof(iint));
  iint *recvOffsets = (iint*) calloc(size+1, sieof(iint));
  for(iint r=0;r<size;++r){
    sendCounts[r] *= sizeof(nonZero_t);
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
  }

  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(sendCounts, 1, MPI_IINT, recvCounts, 1, MPI_IINT, MPI_COMM_WORLD);
  
  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
    sendOffsets[r+1] = sendOffsets[r] + sendCounts[r];
    recvOffsets[r+1] = recvOffsets[r] + recvCounts[r];
    recvNtotal += recvCounts[r]/sizeof(nonZero_t);
  }

  nonZeros_t *recvNonZeros =  (nonZeros_t*) calloc(recvNtotal, sizeof(nonZeros_t));
  MPI_Alltoallv(sendNonZeros, sendCounts, sendOffsets, MPI_CHAR,
		recvNonZeros, recvCounts, recvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);  

  // ------------------------------------------------------------------------------------
  // 3. Sort received non-zeros based on row then column indices
    
  // sort by row then column
  qsort(recvNonZeros, recvNtotal, sizeof(nonZero_t), compareRowCol);

  // compress duplicate entries
  iint cnt = 0;
  for(iint n=1;n<recvNtotal;++n){
    if(!compareRowCol(recvNonZeros+n, recvNonZeros+n-1)){ // match
      recvNonZeros[cnt].val += recvNonZeros[n].val;
    }
    else{
      ++cnt;
      recvNonZeros[cnt] = recvNonZeros[n];
    }
  }
  recvNtotal = cnt+1;

  // ------------------------------------------------------------------------------------
  // 4. Find number of nodes owned by all ranks
  
  iint *allNowned = (iint*) calloc(size, sizeof(iint));
  iint *allStarts = (iint*) calloc(size+1, sizeof(iint));
  MPI_Allgather(&recvNtotal, 1, MPI_IINT, allNowned, 1, MPI_IINT, MPI_COMM_WORLD);

  for(iint r=0;r<size;++r){
    allStarts[r+1] = allStarts[r] + allNowned[r];
  }

  // ------------------------------------------------------------------------------------
  // 5. Set up AMG
  
  precon->amg = amgSetup(allStarts, recvNonZeros);

  free(recvNonZeros);
 
}
