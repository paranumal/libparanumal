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
       iint   nullSpace,            // null space flag
       hgs_t *hgs,                  // gs op for problem assembly (to be removed in future)
       const char* options)      
{
  iint size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iint numLocalRows = globalRowStarts[rank+1]-globalRowStarts[rank];
  iint globalOffset = globalRowStarts[rank];
  iint numGlobalRows = globalRowStarts[size];

  //count number of local, and non-local non-zeros
  iint diag_nnz=0; 
  iint offd_nnz=0; 
  for (iint n=0;n<nnz;n++) {
    if ((Aj[n] < globalOffset) || (Aj[n]>globalOffset+numLocalRows-1)) offd_nnz++;
    else diag_nnz++;
  }

  iint   *diagRowStarts, *diagAi, *diagAj;
  iint   *offdRowStarts, *offdAi, *offdAj;
  dfloat *diagAvals, *offdAvals;

  if (diag_nnz) {
    diagRowStarts = (iint *)   calloc(numLocalRows+1, sizeof(iint));  
    diagAi        = (iint *)   calloc(diag_nnz, sizeof(iint));
    diagAj        = (iint *)   calloc(diag_nnz, sizeof(iint));
    diagAvals     = (dfloat *) calloc(diag_nnz, sizeof(dfloat));
  }
  if (offd_nnz) {
    offdRowStarts = (iint *)   calloc(numLocalRows+1, sizeof(iint));  
    offdAi        = (iint *)   calloc(offd_nnz, sizeof(iint));
    offdAj        = (iint *)   calloc(offd_nnz, sizeof(iint));
    offdAvals     = (dfloat *) calloc(offd_nnz, sizeof(dfloat));
  }

  //split into local and non-local COO matrices
  diag_nnz =0;
  offd_nnz =0;
  for (iint n=0;n<nnz;n++) {
    if ((Aj[n] < globalOffset) || (Aj[n]>globalOffset+numLocalRows-1)) {
      offdAi[offd_nnz] = Ai[n] - globalOffset; //local index
      offdAj[offd_nnz] = Aj[n];                //global index
      offdAvals[offd_nnz] = Avals[n];
      offd_nnz++;
    } else {
      diagAi[diag_nnz] = Ai[n] - globalOffset; //local index
      diagAj[diag_nnz] = Aj[n] - globalOffset; //local index
      diagAvals[diag_nnz] = Avals[n];
      diag_nnz++;
    }
  }

  // Convert to csr storage, assumes orginal matrix was presorted by rows
  iint cnt =0; //row start counter
  for(iint n=0;n<diag_nnz;++n)
    if(cnt==0 || (diagAi[n]!=diagAi[n-1])) 
      diagRowStarts[cnt++] = n;
  diagRowStarts[cnt] = diag_nnz;
  
  cnt =0; //row start counter
  for(iint n=0;n<offd_nnz;++n)
    if(cnt==0 || (offdAi[n]!=offdAi[n-1])) 
      offdRowStarts[cnt++] = n;
  offdRowStarts[cnt] = offd_nnz;

  csr *A = newCSR(numLocalRows, diag_nnz, diagRowStarts, diagAj, diagAvals,
                                offd_nnz, offdRowStarts, offdAj, offdAvals);
  free(diagRowStarts); free(diagAi); free(diagAj); free(diagAvals);
  free(offdRowStarts); free(offdAi); free(offdAj); free(offdAvals);

  //populate null space vector
  dfloat *nullA = (dfloat *) calloc(numLocalRows, sizeof(dfloat));
  for (iint i=0;i<numLocalRows;i++) nullA[i] = 1;
  
  parAlmond_t *parAlmond = agmgSetup(A, nullA, globalRowStarts, options);
  
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

  //sync_setup_on_device(parAlmond, mesh->device);
  return (void *) parAlmond;
}

//TODO code this
int parAlmondFree(void* A) {
  return 0;
}