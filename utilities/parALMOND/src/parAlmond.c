#include "parAlmond.h"

void parAlmondPrecon(occa::memory o_x, void *A, occa::memory o_rhs) {

  parAlmond_t *parAlmond = (parAlmond_t*) A;

  //gather the global problem
  //if the rhs has already been gather scattered, weight the gathered rhs
  if(strstr(parAlmond->options,"CONTINUOUS")||strstr(parAlmond->options,"PROJECT")
      ||strstr(parAlmond->options,"PRECONC0")) {
    meshParallelGather(parAlmond->mesh, parAlmond->hgs, o_rhs, parAlmond->levels[0]->o_rhs);
    dotStar(parAlmond, parAlmond->hgs->Ngather,
            parAlmond->hgs->o_invDegree, parAlmond->levels[0]->o_rhs);
  } else {
    parAlmond->levels[0]->o_rhs.copyFrom(o_rhs);
  }

  device_kcycle(parAlmond, 0);
  //device_vcycle(parAlmond, 0);
  //device_pcg(parAlmond,1000,1e-8);

  //host versions
  //parAlmond->levels[0]->o_rhs.copyTo(parAlmond->levels[0]->rhs);
  //kcycle(parAlmond,0);
  //vcycle(parAlmond,0);
  //pcg(parAlmond,1000,1e-8);
  //parAlmond->levels[0]->o_x.copyFrom(parAlmond->levels[0]->x);

  //scatter the result
  if(strstr(parAlmond->options,"CONTINUOUS")||strstr(parAlmond->options,"PROJECT")
      ||strstr(parAlmond->options,"PRECONC0")) {
    meshParallelScatter(parAlmond->mesh, parAlmond->hgs, parAlmond->levels[0]->o_x, o_x);
  } else {
    parAlmond->levels[0]->o_x.copyTo(o_x);
  }
}

void *parAlmondInit(mesh_t *mesh, const char* options) {

  parAlmond_t *parAlmond = (parAlmond_t *) calloc(1,sizeof(parAlmond_t));

  parAlmond->mesh = mesh; //TODO parALmond doesnt need mesh, except for GS kernels.
  parAlmond->device = mesh->device;
  parAlmond->options = options;

  parAlmond->levels = (agmgLevel **) calloc(MAX_LEVELS,sizeof(agmgLevel *));
  parAlmond->numLevels = 0;
  parAlmond->ktype = PCG;

  buildAlmondKernels(parAlmond);

  //buffer for innerproducts in kcycle
  parAlmond->o_rho  = mesh->device.malloc(3*sizeof(dfloat));

  return (void *) parAlmond;
}

void parAlmondAddDeviceLevel(void *Almond, iint lev, iint Nrows, iint Ncols,
        void **AxArgs,        void (*Ax)(void **args, occa::memory o_x, occa::memory o_Ax),
        void **coarsenArgs,   void (*coarsen)(void **args, occa::memory o_x, occa::memory o_Rx),
        void **prolongateArgs,void (*prolongate)(void **args, occa::memory o_x, occa::memory o_Px),
        void **smootherArgs,  void (*smooth)(void **args, occa::memory o_r, occa::memory o_x, bool x_is_zero)) {

  parAlmond_t *parAlmond = (parAlmond_t *) Almond;
  agmgLevel **levels = parAlmond->levels;

  if (lev > parAlmond->numLevels-1)
    parAlmond->numLevels = lev+1;

  if (!levels[lev])
    levels[lev] = (agmgLevel *) calloc(1,sizeof(agmgLevel));

  levels[lev]->Nrows = Nrows;
  levels[lev]->Ncols = Ncols;

  levels[lev]->AxArgs = AxArgs;
  levels[lev]->device_Ax = Ax;

  levels[lev]->smootherArgs = smootherArgs;
  levels[lev]->device_smooth = smooth;

  if (lev > 0) { //if adding the first level, ignore the coarsen/prologagtion
    levels[lev]->coarsenArgs = coarsenArgs;
    levels[lev]->device_coarsen = coarsen;

    levels[lev]->prolongateArgs = prolongateArgs;
    levels[lev]->device_prolongate = prolongate;
  }

  parAlmond->numLevels++;
}

void parAlmondAgmgSetup(void *Almond,
       iint* globalRowStarts,       //global partition
       iint  nnz,                   //--
       iint* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
       iint* Aj,                    //--
       dfloat* Avals,               //--
       hgs_t *hgs,                  // gs op for problem assembly (to be removed in future?)
       const char* options){

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  parAlmond_t *parAlmond = (parAlmond_t *) Almond;

  iint numLocalRows = globalRowStarts[rank+1]-globalRowStarts[rank];

  csr *A = newCSRfromCOO(numLocalRows,globalRowStarts,nnz, Ai, Aj, Avals);

  //populate null space vector
  dfloat *nullA = (dfloat *) calloc(numLocalRows, sizeof(dfloat));
  for (iint i=0;i<numLocalRows;i++) nullA[i] = 1;

  agmgSetup(parAlmond,A, nullA, globalRowStarts, options);

  parAlmond->hgs = hgs;

  mesh_t *mesh = parAlmond->mesh;

  //buffers for matrix-free action
  parAlmond->o_x = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat));
  parAlmond->o_Ax = mesh->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*sizeof(dfloat));

  if (strstr(options, "VERBOSE"))
    parAlmondReport(parAlmond);
}




//TODO code this
int parAlmondFree(void* A) {
  return 0;
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
