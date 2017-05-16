#include "ellipticQuad2D.h"

typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

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


void ellipticCoarsePreconditionerSetupQuad2D(mesh_t *mesh, precon_t *precon, dfloat lambda, const char *options){

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // ------------------------------------------------------------------------------------
  // 1. Create a contiguous numbering system, starting from the element-vertex connectivity
  iint Nnum = mesh->Nverts*(mesh->Nelements);
  
  iint *globalNumbering = (iint*) calloc(Nnum, sizeof(iint));

  iint *globalOwners = (iint*) calloc(Nnum, sizeof(iint));
  iint *globalStarts = (iint*) calloc(size+1, sizeof(iint));
  
  // use original vertex numbering
  memcpy(globalNumbering, mesh->EToV, mesh->Nelements*mesh->Nverts*sizeof(iint));
  
  // squeeze numbering
  meshParallelConsecutiveGlobalNumbering(Nnum, globalNumbering, globalOwners, globalStarts);
  
  //use the ordering to define a gather+scatter for assembly
  hgs_t *hgs = meshParallelGatherSetup(mesh, Nnum, globalNumbering, globalOwners);


  // temporary
  precon->o_ztmp = mesh->device.malloc(mesh->Np*(mesh->Nelements+mesh->totalHaloPairs)*sizeof(dfloat));
  
  // ------------------------------------------------------------------------------------
  // 2. Build coarse grid element basis functions

  dfloat *V1  = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vr1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));
  dfloat *Vs1 = (dfloat*) calloc(mesh->Np*mesh->Nverts, sizeof(dfloat));

  for(iint j=0;j<mesh->Nq;++j){
    for(iint i=0;i<mesh->Nq;++i){
      iint n = i+j*mesh->Nq;

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
  precon->o_V1  = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), V1);
  precon->o_Vr1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vr1);
  precon->o_Vs1 = mesh->device.malloc(mesh->Nverts*mesh->Np*sizeof(dfloat), Vs1);

  // ------------------------------------------------------------------------------------
  // 3. Build non-zeros of stiffness matrix (unassembled)
  iint nnz = mesh->Nverts*mesh->Nverts*(mesh->Nelements+mesh->totalHaloPairs);
  iint   *rowsA = (iint*) calloc(nnz, sizeof(iint));
  iint   *colsA = (iint*) calloc(nnz, sizeof(iint));
  dfloat *valsA = (dfloat*) calloc(nnz, sizeof(dfloat));

  nonZero_t *sendNonZeros = (nonZero_t*) calloc(nnz, sizeof(nonZero_t));
  iint *AsendCounts  = (iint*) calloc(size, sizeof(iint));
  iint *ArecvCounts  = (iint*) calloc(size, sizeof(iint));
  iint *AsendOffsets = (iint*) calloc(size+1, sizeof(iint));
  iint *ArecvOffsets = (iint*) calloc(size+1, sizeof(iint));
  
  iint cnt = 0;

  printf("Building coarse matrix system\n");
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      for(iint m=0;m<mesh->Nverts;++m){
      	dfloat Snm = 0;
       
      	// use GLL nodes for integration
      	// (since Jacobian is high order tensor-product polynomial)
      	for(iint j=0;j<mesh->Nq;++j){
      	  for(iint i=0;i<mesh->Nq;++i){
      	    iint id = i+j*mesh->Nq;
            
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
  for(iint n=0;n<cnt;++n)
    AsendCounts[sendNonZeros[n].ownerRank] += sizeof(nonZero_t);

  // sort by row ordering
  qsort(sendNonZeros, cnt, sizeof(nonZero_t), parallelCompareRowColumn);
  
  // find how many nodes to expect (should use sparse version)
  MPI_Alltoall(AsendCounts, 1, MPI_IINT, ArecvCounts, 1, MPI_IINT, MPI_COMM_WORLD);

  // find send and recv offsets for gather
  iint recvNtotal = 0;
  for(iint r=0;r<size;++r){
    AsendOffsets[r+1] = AsendOffsets[r] + AsendCounts[r];
    ArecvOffsets[r+1] = ArecvOffsets[r] + ArecvCounts[r];
    recvNtotal += ArecvCounts[r]/sizeof(nonZero_t);
  }

  nonZero_t *recvNonZeros = (nonZero_t*) calloc(recvNtotal, sizeof(nonZero_t));
  
  // determine number to receive
  MPI_Alltoallv(sendNonZeros, AsendCounts, AsendOffsets, MPI_CHAR,
		recvNonZeros, ArecvCounts, ArecvOffsets, MPI_CHAR,
		MPI_COMM_WORLD);

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort(recvNonZeros, recvNtotal, sizeof(nonZero_t), parallelCompareRowColumn);

  // compress duplicates
  cnt = 0;
  for(iint n=1;n<recvNtotal;++n){
    if(recvNonZeros[n].row == recvNonZeros[cnt].row &&
       recvNonZeros[n].col == recvNonZeros[cnt].col){
      recvNonZeros[cnt].val += recvNonZeros[n].val;
    }
    else{
      ++cnt;
      recvNonZeros[cnt] = recvNonZeros[n];
    }
  }
  recvNtotal = cnt+1;

  iint *recvRows = (iint *) calloc(recvNtotal,sizeof(iint));
  iint *recvCols = (iint *) calloc(recvNtotal,sizeof(iint));
  dfloat *recvVals = (dfloat *) calloc(recvNtotal,sizeof(dfloat));
  
  for (iint n=0;n<recvNtotal;n++) {
    recvRows[n] = recvNonZeros[n].row;
    recvCols[n] = recvNonZeros[n].col;
    recvVals[n] = recvNonZeros[n].val;
  }

  //collect global assembled matrix
  iint *globalnnz       = (iint *) calloc(size  ,sizeof(iint));
  iint *globalnnzOffset = (iint *) calloc(size+1,sizeof(iint));
  MPI_Allgather(&recvNtotal, 1, MPI_IINT, 
                globalnnz, 1, MPI_IINT, MPI_COMM_WORLD);
  globalnnzOffset[0] = 0;
  for (iint n=0;n<size;n++)
    globalnnzOffset[n+1] = globalnnzOffset[n]+globalnnz[n];

  iint globalnnzTotal = globalnnzOffset[size];

  iint *globalRecvCounts  = (iint *) calloc(size,sizeof(iint));
  iint *globalRecvOffsets = (iint *) calloc(size,sizeof(iint));
  for (iint n=0;n<size;n++){
    globalRecvCounts[n] = globalnnz[n]*sizeof(nonZero_t);
    globalRecvOffsets[n] = globalnnzOffset[n]*sizeof(nonZero_t);
  }
  nonZero_t *globalNonZero = (nonZero_t*) calloc(globalnnzTotal, sizeof(nonZero_t));

  MPI_Allgatherv(recvNonZeros, recvNtotal*sizeof(nonZero_t), MPI_CHAR, 
                globalNonZero, globalRecvCounts, globalRecvOffsets, MPI_CHAR, MPI_COMM_WORLD);
  

  iint *globalIndex = (iint *) calloc(globalnnzTotal, sizeof(iint));
  iint *globalRows = (iint *) calloc(globalnnzTotal, sizeof(iint));
  iint *globalCols = (iint *) calloc(globalnnzTotal, sizeof(iint));
  dfloat *globalVals = (dfloat*) calloc(globalnnzTotal,sizeof(dfloat));

  for (iint n=0;n<globalnnzTotal;n++) {
    globalRows[n] = globalNonZero[n].row;
    globalCols[n] = globalNonZero[n].col;
    globalVals[n] = globalNonZero[n].val;
  }
  
  if(strstr(options, "XXT")){

    precon->xxt = xxtSetup(Nnum,
			   globalNumbering,
			   nnz,
			   rowsA,
			   colsA,
			   valsA,
			   0,
			   iintString,
			   dfloatString); // 0 if no null space
    
  }

  if(strstr(options, "ALMOND")){
 
    precon->parAlmond = parAlmondSetup(mesh,
         Nnum, 
				 globalStarts,
				 recvNtotal,      
				 recvRows,        
				 recvCols,       
				 recvVals,   
				 0, // 0 if no null space
         hgs,
         options); //rhs will be passed gather-scattered
    
  }

  iint NnumWHalo = Nnum + mesh->Nverts*mesh->totalHaloPairs;
  precon->o_r1 = mesh->device.malloc(NnumWHalo*sizeof(dfloat));
  precon->o_z1 = mesh->device.malloc(NnumWHalo*sizeof(dfloat));
  precon->r1 = (dfloat*) malloc(NnumWHalo*sizeof(dfloat));
  precon->z1 = (dfloat*) malloc(NnumWHalo*sizeof(dfloat));
  

  free(AsendCounts);
  free(ArecvCounts);
  free(AsendOffsets);
  free(ArecvOffsets);

}
