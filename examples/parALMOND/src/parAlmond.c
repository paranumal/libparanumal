#include "agmg.h"

void parAlmondPrecon(parAlmond_t *parAlmond, occa::memory o_x, occa::memory o_rhs) {
  //gather the global problem
  //if the rhs has already been gather scattered, weight the gathered rhs
  if(strstr(parAlmond->options,"GATHER")) {
    meshParallelGather(parAlmond->mesh, parAlmond->hgs, o_rhs, parAlmond->levels[0]->o_rhs);
    dotStar(parAlmond, parAlmond->hgs->Ngather,
            parAlmond->hgs->o_invDegree, parAlmond->levels[0]->o_rhs);
  } else {
    parAlmond->levels[0]->o_rhs.copyFrom(o_rhs);
  }

  if (strstr(parAlmond->options,"HOST")) {
    //host versions
    parAlmond->levels[0]->o_rhs.copyTo(parAlmond->levels[0]->rhs);
    if(strstr(parAlmond->options,"KCYCLE")) {
      kcycle(parAlmond, 0);
    } else if(strstr(parAlmond->options,"VCYCLE")) {
      vcycle(parAlmond, 0);
    } else if(strstr(parAlmond->options,"EXACT")) {
      pcg(parAlmond,1000,1e-8);
    }
    parAlmond->levels[0]->o_x.copyFrom(parAlmond->levels[0]->x);
  } else {
    if(strstr(parAlmond->options,"KCYCLE")) {
      device_kcycle(parAlmond, 0);
    } else if(strstr(parAlmond->options,"VCYCLE")) {
      device_vcycle(parAlmond, 0);
    } else if(strstr(parAlmond->options,"EXACT")) {
      device_pcg(parAlmond,1000,1e-8);
    }
  }

  //scatter the result
  if(strstr(parAlmond->options,"GATHER")) {
    meshParallelScatter(parAlmond->mesh, parAlmond->hgs, parAlmond->levels[0]->o_x, o_x);
  } else {
    parAlmond->levels[0]->o_x.copyTo(o_x);
  }
}

parAlmond_t *parAlmondInit(mesh_t *mesh, const char* options) {

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

  return parAlmond;
}

void parAlmondAgmgSetup(parAlmond_t *parAlmond,
       iint* globalRowStarts,       //global partition
       unsigned long long  nnz,                   //--
       iint* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
       iint* Aj,                    //--
       dfloat* Avals,               //--
       hgs_t *hgs){                  // gs op for problem assembly (to be removed in future?)

  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint numLocalRows = globalRowStarts[rank+1]-globalRowStarts[rank];

  csr *A = newCSRfromCOO(numLocalRows,globalRowStarts,nnz, Ai, Aj, Avals);

  //populate null space vector
  dfloat *nullA = (dfloat *) calloc(numLocalRows, sizeof(dfloat));
  for (iint i=0;i<numLocalRows;i++) nullA[i] = 1;

  agmgSetup(parAlmond, A, nullA, globalRowStarts, parAlmond->options);

  parAlmond->hgs = hgs;

  mesh_t *mesh = parAlmond->mesh;

  if (strstr(parAlmond->options, "VERBOSE"))
    parAlmondReport(parAlmond);
}

//TODO code this
int parAlmondFree(void* A) {
  return 0;
}

