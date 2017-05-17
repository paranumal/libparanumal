
#include "parAlmond.h"


void * parAlmondSetup(mesh_t *mesh, 
       iint  Nnum,
       iint* rowStarts, 
       iint  nnz, 
       iint* Ai,
       iint* Aj,
       dfloat* Avals,
       iint   nullSpace,
       hgs_t *hgs,
       const char* options)
{
  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint numLocalRows = rowStarts[rank+1]-rowStarts[rank];
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
  
  csr *A = newCSR(numGlobalRows, numGlobalRows, nnz, 
                        vRowStarts, vAj, vAvals);
  free(vRowStarts);
  free(vAj);
  free(vAvals);

  dfloat *nullA = (dfloat *) calloc(numGlobalRows, sizeof(dfloat));
  for (iint i=0;i<numGlobalRows;i++) nullA[i] = 1;
  
  parAlmond_t *parAlmond = agmgSetup(A, nullA, rowStarts, options);
  
  
  freeCSR(A);
  free(nullA);

  parAlmond->mesh = mesh;
  parAlmond->hgs = hgs;
  parAlmond->options = options;
  parAlmond->ktype = PCG;

  //buffers for matrix-free action
  dfloat *dummy = (dfloat *) malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat));
  parAlmond->o_x = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat),dummy);
  parAlmond->o_Ax = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat),dummy);
  free(dummy);

sync_setup_on_device(parAlmond, mesh->device);
  return (void *) parAlmond;
}

void parAlmondPrecon(occa::memory o_x, void *A, occa::memory o_rhs) {

  parAlmond_t *parAlmond = (parAlmond_t*) A;

  //gather the global problem
  meshParallelGather(parAlmond->mesh, parAlmond->hgs, o_rhs, parAlmond->levels[0]->o_rhs);

  //if the rhs has already been gather scattered, weight the gathered rhs
  if(strstr(parAlmond->options,"CONTINUOUS")||strstr(parAlmond->options,"PROJECT")) {
    parAlmond->mesh->dotMultiplyKernel(parAlmond->hgs->Ngather,parAlmond->hgs->o_invDegree,
                          parAlmond->levels[0]->o_rhs, parAlmond->levels[0]->o_rhs);
  }

  //parAlmond->levels[0]->o_rhs.copyTo(parAlmond->levels[0]->rhs);
  ////kcycle(parAlmond,0);
  //vcycle(parAlmond,0);
  //parAlmond->levels[0]->o_x.copyFrom(parAlmond->levels[0]->x);

  device_kcycle(parAlmond, 0);
  //device_vcycle(parAlmond, 0);

  //scatter the result
  meshParallelScatter(parAlmond->mesh, parAlmond->hgs, parAlmond->levels[0]->o_x, o_x);
}

//TODO code this
int parAlmondFree(void* A) {
  return 0;
}

void parAlmondMatrixFreeAX(parAlmond_t *parAlmond, occa::memory &o_x, occa::memory &o_Ax){

  mesh_t* mesh = parAlmond->mesh;

  //scatter x 
  meshParallelScatter(mesh, parAlmond->hgs, o_x, parAlmond->o_x);

  occaTimerTic(mesh->device,"MatFreeAxKernel");
  parAlmond->MatFreeAx(parAlmond->MatFreeArgs,parAlmond->o_x,parAlmond->o_Ax,parAlmond->options);
  occaTimerToc(mesh->device,"MatFreeAxKernel");

  //gather the result back to the global problem
  meshParallelGather(mesh, parAlmond->hgs, parAlmond->o_Ax, o_Ax);
}

void parAlmondSetMatFreeAX(void* A, void (*MatFreeAx)(void **args, occa::memory o_q, occa::memory o_Aq,const char* options),
                        void **args) {
  parAlmond_t *parAlmond = (parAlmond_t*) A;

  parAlmond->MatFreeAx = MatFreeAx;
  parAlmond->MatFreeArgs = args;

}

