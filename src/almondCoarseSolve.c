
#include "parAlmond.h"

typedef struct {
  csr *A;
  almond_t *almond;

  hgs_t *hgs;
  iint preGathered;

  iint numLocalRows;
  iint Nnum;

  occa::memory o_rhs, o_x;

} parAlmond_t;

void * almondSetup(mesh_t *mesh, 
       iint  Nnum,
       iint* rowStarts, 
       iint  nnz, 
       iint* Ai,
       iint* Aj,
       dfloat* Avals,
       iint   nullSpace,
       hgs_t *hgs,
       iint preGathered) //1 if the rhs is gather-scattered before calling solve
{

  parAlmond_t *parAlmond = (parAlmond_t*) calloc(1, sizeof(parAlmond_t));

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  parAlmond->hgs = hgs;
  parAlmond->preGathered = preGathered;

  parAlmond->numLocalRows = rowStarts[rank+1]-rowStarts[rank];
  
  iint globalOffset = rowStarts[rank];

  iint numGlobalRows = rowStarts[size];

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

  dfloat *nullA = (dfloat *) calloc(numGlobalRows, sizeof(dfloat));
  for (iint i=0;i<numGlobalRows;i++) nullA[i] = 1;
  
  parAlmond->almond = setup(parAlmond->A, nullA, rowStarts);
  sync_setup_on_device(parAlmond->almond, mesh->device);
  
  free(vRowStarts);
  free(vAj);
  free(vAvals);
  free(nullA);

  parAlmond->almond->ktype = PCG;
  
  dfloat *dummy = (dfloat *) calloc(parAlmond->numLocalRows,sizeof(dfloat));
  parAlmond->o_rhs = mesh->device.malloc(parAlmond->numLocalRows*sizeof(dfloat), dummy);
  parAlmond->o_x   = mesh->device.malloc(parAlmond->numLocalRows*sizeof(dfloat), dummy);
  free(dummy);

  return (void *) parAlmond;
}

void almondSolve(mesh_t *mesh, occa::memory o_x, void *A, occa::memory o_rhs) {

  parAlmond_t *parAlmond = (parAlmond_t*) A;

  //gather the global problem
  meshParallelGather(mesh, parAlmond->hgs, o_rhs, parAlmond->o_rhs);

  //if the rhs has already been gather scattered, weight the gathered rhs
  if(parAlmond->preGathered)
    mesh->dotMultiplyKernel(parAlmond->hgs->Ngather,parAlmond->hgs->o_invDegree,
                          parAlmond->o_rhs, parAlmond->o_rhs);

  solve(parAlmond->almond, parAlmond->o_rhs, parAlmond->o_x);

  //scatter the result
  meshParallelScatter(mesh, parAlmond->hgs, parAlmond->o_x, o_x);
}

//TODO code this
int almondFree(void* A) {
  return 0;
}

void almondCoarseSolveSetup(void *ALMOND, iint *coarseNp, iint *coarseOffsets, iint **globalNumbering,
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

/*
void almondMatFreeAX(occa::memory &o_q, occa::memory &o_Aq){

  occaTimerTic(mesh->device,"MatFreeAxKernel");
  
  dfloat *sendBuffer = almond->sendBuffer;
  dfloat *recvBuffer = almond->recvBuffer;

  // compute local element operations and store result in o_Aq
  if(strstr(options, "CONTINUOUS")){
    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);

    // parallel gather scatter
    ellipticParallelGatherScatterHex3D(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");
  }
  else{

    ellipticStartHaloExchange3D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_q, recvBuffer);

    iint allNelements = mesh->Nelements+mesh->totalHaloPairs;
    mesh->gradientKernel(allNelements, mesh->o_vgeo, mesh->o_D, o_q, solver->o_grad);

    dfloat tau = 2.f*(mesh->Nq)*(mesh->Nq+2)/3.;
    mesh->ipdgKernel(mesh->Nelements,
         mesh->o_vmapM,
         mesh->o_vmapP,
         lambda,
         tau,
         mesh->o_vgeo,
         mesh->o_sgeo,
         mesh->o_D,
         solver->o_grad,
         o_Aq);
        
  }

  if(strstr(options, "CONTINUOUS")||strstr(options, "PROJECT"))
    // parallel gather scatter
    ellipticParallelGatherScatterHex3D(mesh, solver->ogs, o_Aq, o_Aq, dfloatString, "add");
 
  occaTimerToc(mesh->device,"AxKernel");
}

void almondSetMatFreeAX(void* ALMOND, ,
                          void *ACoarse, iint coarseTotal,iint coarseOffset) {

}
*/
