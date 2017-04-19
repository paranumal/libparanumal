
#include "parAlmond.h"


typedef struct {
  csr *A;
  almond_t *almond;
  dfloat *rhs;
  dfloat *x;
  dfloat *nullA;

  iint numLocalRows;
  iint Nnum;
  iint recvNnum;
  iint nnz;

  iint *sendSortId, *globalSortId, *compressId;
  iint *sendCounts, *sendOffsets,  *recvCounts, *recvOffsets;

  dfloat *xUnassembled;
  dfloat *rhsUnassembled;

  dfloat *xSort;
  dfloat *rhsSort;

  occa::memory o_rhs, o_x;

} parAlmond_t;

void * almondSetup(occa::device device,
       iint  Nnum,
		   iint* rowStarts, 
		   iint  nnz, 
		   iint* Ai,
		   iint* Aj,
		   dfloat* Avals,
		   iint    *sendSortId, 
		   iint    *globalSortId, 
		   iint    *compressId,
		   iint    *sendCounts, 
		   iint    *sendOffsets, 
		   iint    *recvCounts, 
		   iint    *recvOffsets,
		   iint   nullSpace) {

  parAlmond_t *parAlmond = (parAlmond_t*) calloc(1, sizeof(parAlmond_t));

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  parAlmond->numLocalRows = rowStarts[rank+1]-rowStarts[rank];
  iint globalOffset = rowStarts[rank];

  iint *vRowStarts = (iint *) calloc(parAlmond->numLocalRows+1, sizeof(iint));

  // assumes presorted
  iint cnt = 0; //nnz counter
  iint cnt2 =0; //row start counter
  for(iint n=0;n<nnz;++n){
    if(  (Ai[n] >= (parAlmond->numLocalRows + globalOffset)) || (Ai[n] < globalOffset))
      printf("errant nonzero %d,%d,%g, rank %d \n", Ai[n], Aj[n], Avals[n], rank);

    if(cnt2==0 || (Ai[n]!=Ai[n-1])) vRowStarts[cnt2++] = cnt;
      
    if (Aj[n] >= globalOffset && Aj[n] < parAlmond->numLocalRows + globalOffset) cnt++;
  }

  vRowStarts[cnt2] = cnt;

  iint *vAj = (iint *) calloc(cnt, sizeof(iint));
  dfloat *vAvals = (dfloat *) calloc(cnt, sizeof(dfloat));

  cnt = 0;
  for(iint n=0;n<nnz;++n){
    if (Aj[n] >= globalOffset && Aj[n] < parAlmond->numLocalRows + globalOffset) {
      vAj[cnt] = Aj[n]-globalOffset;
      vAvals[cnt++] = Avals[n];  
    }
  }

  parAlmond->Nnum = Nnum;
  //parAlmond->numLocalRows = numLocalRows;
  parAlmond->recvNnum = compressId[parAlmond->numLocalRows];

  parAlmond->sendSortId = (iint*) calloc(Nnum,sizeof(iint));
  for (iint n=0;n<Nnum;n++)  parAlmond->sendSortId[n] = sendSortId[n];

  parAlmond->sendCounts = (iint*) calloc(size,sizeof(iint));
  parAlmond->recvCounts = (iint*) calloc(size,sizeof(iint));
  for (iint n=0;n<size;n++){
    parAlmond->sendCounts[n] = sendCounts[n]*sizeof(dfloat);
    parAlmond->recvCounts[n] = recvCounts[n]*sizeof(dfloat);
  }

  parAlmond->sendOffsets = (iint*) calloc(size+1,sizeof(iint));
  parAlmond->recvOffsets = (iint*) calloc(size+1,sizeof(iint));
  for (iint n=0;n<size+1;n++){
    parAlmond->sendOffsets[n] = sendOffsets[n]*sizeof(dfloat);
    parAlmond->recvOffsets[n] = recvOffsets[n]*sizeof(dfloat);
  }

  parAlmond->globalSortId = (iint*) calloc(parAlmond->recvNnum,sizeof(iint));
  for (iint n=0;n<parAlmond->recvNnum;n++) parAlmond->globalSortId[n] = globalSortId[n];

  parAlmond->compressId  = (iint*) calloc(parAlmond->numLocalRows+1,sizeof(iint));
  for (iint n=0;n<parAlmond->numLocalRows+1;n++) parAlmond->compressId[n] = compressId[n];

  iint Nmax =  (Nnum >parAlmond->recvNnum)? Nnum : parAlmond->recvNnum;
  parAlmond->xUnassembled   = (dfloat *) malloc(Nmax*sizeof(dfloat));
  parAlmond->rhsUnassembled = (dfloat *) malloc(Nmax*sizeof(dfloat));
  parAlmond->xSort   = (dfloat *) malloc(Nmax*sizeof(dfloat));
  parAlmond->rhsSort = (dfloat *) malloc(Nmax*sizeof(dfloat));

  parAlmond->rhs = (dfloat *) calloc(parAlmond->numLocalRows, sizeof(dfloat));
  parAlmond->x   = (dfloat *) calloc(parAlmond->numLocalRows, sizeof(dfloat));

  parAlmond->A = newCSR(parAlmond->numLocalRows, parAlmond->numLocalRows, nnz, 
                        vRowStarts, vAj, vAvals);

  parAlmond->nullA = (dfloat *) calloc(parAlmond->numLocalRows, sizeof(dfloat));
  for (iint i=0;i<parAlmond->numLocalRows;i++) parAlmond->nullA[i] = 1;
  
  parAlmond->almond = setup(parAlmond->A, parAlmond->nullA, NULL);
  sync_setup_on_device(parAlmond->almond, device);

  parAlmond->almond->ktype = PCG;
  
  return (void *) parAlmond;
}


void * almondGlobalSetup(occa::device device, 
       iint  Nnum,
       iint* rowStarts, 
       iint  nnz, 
       iint* Ai,
       iint* Aj,
       dfloat* Avals,
       iint    *sendSortId, 
       iint    *globalSortId, 
       iint    *compressId,
       iint    *sendCounts, 
       iint    *sendOffsets, 
       iint    *recvCounts, 
       iint    *recvOffsets,
       iint   nullSpace) {

  parAlmond_t *parAlmond = (parAlmond_t*) calloc(1, sizeof(parAlmond_t));

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  parAlmond->numLocalRows = rowStarts[rank+1]-rowStarts[rank];
  int globalOffset = rowStarts[rank];

  int numGlobalRows = rowStarts[size];

  iint *vRowStarts = (iint *) calloc(numGlobalRows+1, sizeof(iint));  
  iint *vAj        = (iint *) calloc(nnz, sizeof(iint));
  dfloat *vAvals   = (dfloat *) calloc(nnz, sizeof(dfloat));

  // assumes presorted
  iint cnt2 =0; //row start counter
  for(iint n=0;n<nnz;++n) {
    if(cnt2==0 || (Ai[n]!=Ai[n-1])) vRowStarts[cnt2++] = n;      
    vAj[n] = Aj[n];
    vAvals[n] = Avals[n];
  }
  vRowStarts[cnt2] = nnz;
  
  parAlmond->A = newCSR(numGlobalRows, numGlobalRows, nnz, 
                        vRowStarts, vAj, vAvals);

  parAlmond->nullA = (dfloat *) calloc(numGlobalRows, sizeof(dfloat));
  for (iint i=0;i<numGlobalRows;i++) parAlmond->nullA[i] = 1;
  
  parAlmond->almond = setup(parAlmond->A, parAlmond->nullA, rowStarts);
  sync_setup_on_device(parAlmond->almond, device);
  
  parAlmond->almond->ktype = PCG;

  parAlmond->Nnum = Nnum;
  parAlmond->recvNnum = compressId[parAlmond->numLocalRows];

  parAlmond->sendSortId = (iint*) calloc(Nnum,sizeof(iint));
  for (iint n=0;n<Nnum;n++)  parAlmond->sendSortId[n] = sendSortId[n];

  parAlmond->sendCounts = (iint*) calloc(size,sizeof(iint));
  parAlmond->recvCounts = (iint*) calloc(size,sizeof(iint));
  for (iint n=0;n<size;n++){
    parAlmond->sendCounts[n] = sendCounts[n]*sizeof(dfloat);
    parAlmond->recvCounts[n] = recvCounts[n]*sizeof(dfloat);
  }

  parAlmond->sendOffsets = (iint*) calloc(size+1,sizeof(iint));
  parAlmond->recvOffsets = (iint*) calloc(size+1,sizeof(iint));
  for (iint n=0;n<size+1;n++){
    parAlmond->sendOffsets[n] = sendOffsets[n]*sizeof(dfloat);
    parAlmond->recvOffsets[n] = recvOffsets[n]*sizeof(dfloat);
  }

  parAlmond->globalSortId = (iint*) calloc(parAlmond->recvNnum,sizeof(iint));
  for (iint n=0;n<parAlmond->recvNnum;n++) parAlmond->globalSortId[n] = globalSortId[n];

  parAlmond->compressId  = (iint*) calloc(parAlmond->numLocalRows+1,sizeof(iint));
  for (iint n=0;n<parAlmond->numLocalRows+1;n++) parAlmond->compressId[n] = compressId[n];

  iint Nmax =  (Nnum >parAlmond->recvNnum)? Nnum : parAlmond->recvNnum;
  parAlmond->xUnassembled   = (dfloat *) malloc(Nmax*sizeof(dfloat));
  parAlmond->rhsUnassembled = (dfloat *) malloc(Nmax*sizeof(dfloat));
  parAlmond->xSort   = (dfloat *) malloc(Nmax*sizeof(dfloat));
  parAlmond->rhsSort = (dfloat *) malloc(Nmax*sizeof(dfloat));

  parAlmond->rhs = (dfloat *) calloc(parAlmond->numLocalRows, sizeof(dfloat));
  parAlmond->x   = (dfloat *) calloc(parAlmond->numLocalRows, sizeof(dfloat));
  
  parAlmond->o_rhs = device.malloc(parAlmond->numLocalRows*sizeof(dfloat), parAlmond->rhs);
  parAlmond->o_x   = device.malloc(parAlmond->numLocalRows*sizeof(dfloat), parAlmond->x);

  return (void *) parAlmond;
}


iint almondSolve(dfloat* x,
		void* A,
		dfloat* rhs) {
  

  parAlmond_t *parAlmond = (parAlmond_t*) A;

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  for (iint n=0;n<parAlmond->Nnum;n++)
    parAlmond->rhsUnassembled[n] = rhs[n];

  
  //sort by owner
  for (iint n=0;n<parAlmond->Nnum;n++) 
    parAlmond->rhsSort[n] = parAlmond->rhsUnassembled[parAlmond->sendSortId[n]];

  //Scatter nodes to their owners
  MPI_Alltoallv(parAlmond->rhsSort, parAlmond->sendCounts, parAlmond->sendOffsets, MPI_CHAR,
                parAlmond->rhsUnassembled, parAlmond->recvCounts, parAlmond->recvOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort by globalid 
  for (iint n=0;n<parAlmond->recvNnum;n++) 
    parAlmond->rhsSort[n] = parAlmond->rhsUnassembled[parAlmond->globalSortId[n]];

  //gather
  for(iint n=0;n<parAlmond->numLocalRows;++n){
    parAlmond->rhs[n] = 0.;
    parAlmond->x[n] = 0.;
    for (iint id=parAlmond->compressId[n];id<parAlmond->compressId[n+1];id++) 
      parAlmond->rhs[n] += parAlmond->rhsSort[id];
  }

  parAlmond->o_rhs.copyFrom(parAlmond->rhs);

  if(1){
    //solve(parAlmond->almond, parAlmond->rhs, parAlmond->x);
    solve(parAlmond->almond, parAlmond->o_rhs, parAlmond->o_x);
  } else{
    iint maxIt = 40;
    dfloat tol = 1e-2;
    pcg(parAlmond->almond,
        parAlmond->A,
			  parAlmond->rhs,
			  parAlmond->x,
			  maxIt,
			  tol);
  }

  parAlmond->o_x.copyTo(parAlmond->x);

  //scatter
  for(iint n=0;n<parAlmond->numLocalRows;++n){
    for (iint id = parAlmond->compressId[n];id<parAlmond->compressId[n+1];id++) 
      parAlmond->xSort[id] = parAlmond->x[n]; 
  }

  //sort by original rank
  for (iint n=0;n<parAlmond->recvNnum;n++) 
    parAlmond->xUnassembled[parAlmond->globalSortId[n]] = parAlmond->xSort[n];

  //Scatter nodes back to their original rank
  MPI_Alltoallv(parAlmond->xUnassembled, parAlmond->recvCounts, parAlmond->recvOffsets, MPI_CHAR,
                parAlmond->xSort, parAlmond->sendCounts, parAlmond->sendOffsets, MPI_CHAR,
                MPI_COMM_WORLD);

  //sort back to original ordering
  for (iint n=0;n<parAlmond->Nnum;n++) 
    parAlmond->xUnassembled[parAlmond->sendSortId[n]] = parAlmond->xSort[n];


  for(iint i=0;i<parAlmond->Nnum;++i) x[i] = parAlmond->xUnassembled[i];

  return 0;
}

//TODO code this
int almondFree(void* A) {
  return 0;
}


void almondProlongateCoarseProblem(void *ALMOND, iint *coarseNp, iint *coarseOffsets, dfloat **B) { 

  parAlmond_t *parAlmond = (parAlmond_t*) ALMOND;

  iint rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );

  iint numLevels = parAlmond->almond->numLevels;
  
  //number of coarse basis vectors
  iint localCoarseNp = parAlmond->almond->levels[numLevels-1]->A->Nrows;

  MPI_Allgather(&localCoarseNp, 1, MPI_IINT, coarseNp, 1, MPI_IINT, MPI_COMM_WORLD);

  for (iint n=0;n<size;n++) coarseOffsets[n+1] = coarseOffsets[n] + coarseNp[n];

  iint coarseTotal = coarseOffsets[size];

  *B = (dfloat *) calloc(localCoarseNp*parAlmond->Nnum, sizeof(dfloat));

  //allocate storage for prolanged basis vector at all levels
  dfloat **b = (dfloat **) calloc(numLevels,sizeof(dfloat *));
  for (iint n=0; n<numLevels;n++)
    b[n] = (dfloat *) calloc(parAlmond->almond->levels[n]->Ncols, sizeof(dfloat));

  for (iint k=0;k<localCoarseNp;k++) {
    //fill coarsest vector
    for (iint n=0;n<localCoarseNp;n++) b[numLevels-1][n] = 0.0;
    
    b[numLevels-1][k] = 1.0;

    //prolongate
    for (iint m=numLevels-1;m>0;m--) 
      interpolate(parAlmond->almond->levels[m-1], b[m], b[m-1]);
    
    //scatter
    for(iint n=0;n<parAlmond->numLocalRows;++n){
      for (iint id = parAlmond->compressId[n];id<parAlmond->compressId[n+1];id++) 
        parAlmond->xSort[id] = b[0][n]; 
    }

    //sort by original rank
    for (iint n=0;n<parAlmond->recvNnum;n++) 
      parAlmond->xUnassembled[parAlmond->globalSortId[n]] = parAlmond->xSort[n];

    //Fake scatter nodes back to their original rank
    iint recvCount = parAlmond->recvCounts[rank]/sizeof(dfloat);
    iint recvOffset = parAlmond->recvOffsets[rank]/sizeof(dfloat);
    iint sendOffset = parAlmond->sendOffsets[rank]/sizeof(dfloat);
    for (iint n=0; n<parAlmond->Nnum; n++) parAlmond->xSort[n] = 0.;
    for (iint n=0; n<recvCount;n++) parAlmond->xSort[n+sendOffset] = parAlmond->xUnassembled[n+recvOffset];

    //sort back to original ordering
    for (iint n=0;n<parAlmond->Nnum;n++) 
      parAlmond->xUnassembled[parAlmond->sendSortId[n]] = parAlmond->xSort[n];
    
    //save
    for (iint n=0;n<parAlmond->Nnum;n++) 
      (*B)[n+k*parAlmond->Nnum] = parAlmond->xUnassembled[n];
  }

  for (iint n=0; n<numLevels;n++) 
    free(b[n]);
}

void almondGlobalCoarseSetup(void *ALMOND, iint *coarseNp, iint *coarseOffsets, iint **globalNumbering,
                    iint *nnz, iint **rows, iint **cols, dfloat **vals) {

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  parAlmond_t *parAlmond = (parAlmond_t*) ALMOND;

  iint numLevels = parAlmond->almond->numLevels;

  iint Np = parAlmond->almond->levels[numLevels-1]->Nrows;

  MPI_Allgather(&Np, 1, MPI_IINT, coarseNp, 1, MPI_IINT, MPI_COMM_WORLD);

  coarseOffsets[0] = 0;
  for (iint r=0;r<size;r++)
    coarseOffsets[r+1] = coarseOffsets[r] + coarseNp[r];

  *globalNumbering = (iint *) calloc(coarseOffsets[size],sizeof(iint));
  for (iint n=0;n<coarseOffsets[size];n++)
    (*globalNumbering)[n] = n;

  *nnz = parAlmond->almond->levels[numLevels-1]->A->nnz;
  if (*nnz) {
    *rows = (iint *) calloc(*nnz,sizeof(iint));
    *cols = (iint *) calloc(*nnz,sizeof(iint));
    *vals = (dfloat *) calloc(*nnz,sizeof(dfloat));
  }

  for (iint n=0;n<Np;n++) {
    for (iint m=parAlmond->almond->levels[numLevels-1]->A->rowStarts[n];
             m<parAlmond->almond->levels[numLevels-1]->A->rowStarts[n+1];m++) {
      (*rows)[m] = n + coarseOffsets[rank];
      iint col = parAlmond->almond->levels[numLevels-1]->A->cols[m];
      (*cols)[m] = parAlmond->almond->levels[numLevels-1]->A->colMap[col];
      (*vals)[m] = parAlmond->almond->levels[numLevels-1]->A->coefs[m];
    }
  }
}

void almondSetCoarseSolve(void* ALMOND, void (*coarseSolve)(void*,void*,void*),
                          void *ACoarse, iint coarseTotal,iint coarseOffset) {

  parAlmond_t *parAlmond = (parAlmond_t*) ALMOND;

  //set coarse solver pointer
  parAlmond->almond->coarseSolve = coarseSolve;
  parAlmond->almond->ACoarse = ACoarse;
  parAlmond->almond->coarseTotal = coarseTotal;
  parAlmond->almond->coarseOffset = coarseOffset;
}