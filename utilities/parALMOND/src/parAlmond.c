#include "parAlmond.h"

void parAlmondPrecon(occa::memory o_x, void *A, occa::memory o_rhs) {

  parAlmond_t *parAlmond = (parAlmond_t*) A;

  //gather the global problem
  //if the rhs has already been gather scattered, weight the gathered rhs
  if(strstr(parAlmond->options,"CONTINUOUS")||strstr(parAlmond->options,"PROJECT")) {
    meshParallelGather(parAlmond->mesh, parAlmond->hgs, o_rhs, parAlmond->levels[0]->o_rhs);
    dotStar(parAlmond, parAlmond->hgs->Ngather,
            parAlmond->hgs->o_invDegree, parAlmond->levels[0]->o_rhs);
  } else {
    parAlmond->levels[0]->o_rhs.copyFrom(o_rhs);
  }

  device_kcycle(parAlmond, 0);
  //device_vcycle(parAlmond, 0);

  //parAlmond->levels[0]->o_rhs.copyTo(parAlmond->levels[0]->rhs);
  //kcycle(parAlmond, 0);
  //parAlmond->levels[0]->o_x.copyFrom(parAlmond->levels[0]->x);
/*
  iint M = parAlmond->levels[0]->Nrows;
  occa::memory o_x0 = parAlmond->device.malloc(M*sizeof(dfloat));
  occa::memory o_r0 = parAlmond->device.malloc(M*sizeof(dfloat));

  o_r0.copyFrom(parAlmond->levels[0]->o_rhs);
  pcg(parAlmond,o_r0,o_x0,100,1e-7);
  parAlmond->levels[0]->o_x.copyFrom(o_x0);
  o_r0.free();
  o_x0.free();
*/
  //scatter the result
  if(strstr(parAlmond->options,"CONTINUOUS")||strstr(parAlmond->options,"PROJECT")) {
    meshParallelScatter(parAlmond->mesh, parAlmond->hgs, parAlmond->levels[0]->o_x, o_x);
  } else {
    parAlmond->levels[0]->o_x.copyTo(o_x);
  }
}

void parAlmondMatrixFreeAX(parAlmond_t *parAlmond, occa::memory &o_x, occa::memory &o_Ax){

  mesh_t* mesh = parAlmond->mesh;

  if(strstr(parAlmond->options,"CONTINUOUS")||strstr(parAlmond->options,"PROJECT")) {
    //scatter x
    meshParallelScatter(mesh, parAlmond->hgs, o_x, parAlmond->o_x);

    occaTimerTic(mesh->device,"MatFreeAxKernel");
    parAlmond->MatFreeAx(parAlmond->MatFreeArgs,parAlmond->o_x,parAlmond->o_Ax,parAlmond->options);
    occaTimerToc(mesh->device,"MatFreeAxKernel");

    //gather the result back to the global problem
    meshParallelGather(mesh, parAlmond->hgs, parAlmond->o_Ax, o_Ax);
  } else {
    occaTimerTic(mesh->device,"MatFreeAxKernel");
    parAlmond->MatFreeAx(parAlmond->MatFreeArgs,o_x,o_Ax,parAlmond->options);
    occaTimerToc(mesh->device,"MatFreeAxKernel");
  }
}

void parAlmondSetMatFreeAX(void* A, void (*MatFreeAx)(void **args, occa::memory o_q, occa::memory o_Aq,const char* options),
                        void **args) {
  parAlmond_t *parAlmond = (parAlmond_t*) A;

  parAlmond->MatFreeAx = MatFreeAx;
  parAlmond->MatFreeArgs = args;

}

void * parAlmondSetup(mesh_t *mesh, //mesh data
       iint* globalRowStarts,       //global partition
       iint  nnz,                   //--
       iint* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
       iint* Aj,                    //--
       dfloat* Avals,               //--
       hgs_t *hgs,                  // gs op for problem assembly (to be removed in future)
       const char* options)
{
  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint numLocalRows = globalRowStarts[rank+1]-globalRowStarts[rank];

  csr *A = newCSRfromCOO(numLocalRows,globalRowStarts,nnz, Ai, Aj, Avals);

  //populate null space vector
  dfloat *nullA = (dfloat *) calloc(numLocalRows, sizeof(dfloat));
  for (iint i=0;i<numLocalRows;i++) nullA[i] = 1;

  parAlmond_t *parAlmond = agmgSetup(A, nullA, globalRowStarts, options);

  parAlmond->mesh = mesh;
  parAlmond->hgs = hgs;
  parAlmond->options = options;
  parAlmond->ktype = PCG;

  //buffers for matrix-free action
  parAlmond->o_x = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat));
  parAlmond->o_Ax = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat));

  sync_setup_on_device(parAlmond, mesh->device);
  parAlmondReport(parAlmond);
  return (void *) parAlmond;
}

//TODO code this
int parAlmondFree(void* A) {
  return 0;
}