#include "ellipticTri2D.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(solver_t* solver, const char *options){  

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // initialize gather structs
  solver->halo    = (ogs_t*) calloc(1, sizeof(ogs_t));
  solver->nonHalo = (ogs_t*) calloc(1, sizeof(ogs_t));

  mesh_t *mesh = solver->mesh;

  iint Nlocal = mesh->Np*mesh->Nelements;

  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process
  solver->halo->Ngather = 0;
  solver->nonHalo->Ngather = 0;

  iint nHalo = 0;
  iint nNonHalo = 0;

  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (mesh->gatherBaseIds[n] != mesh->gatherBaseIds[n-1]);
    if(mesh->gatherHaloFlags[n]==1){
      solver->halo->Ngather += test;
      ++nHalo;
    }
  }

  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (mesh->gatherBaseIds[n] != mesh->gatherBaseIds[n-1]);
    if(mesh->gatherHaloFlags[n]!=1){
      solver->nonHalo->Ngather += test;
      ++nNonHalo;
    }
  }

  solver->halo->gatherOffsets  = (iint*) calloc(solver->halo->Ngather+1, sizeof(iint));
  solver->halo->gatherLocalIds = (iint*) calloc(nHalo, sizeof(iint));
  solver->halo->gatherBaseIds  = (iint*) calloc(solver->halo->Ngather, sizeof(iint));

  solver->nonHalo->gatherOffsets  = (iint*) calloc(solver->nonHalo->Ngather+1, sizeof(iint));
  solver->nonHalo->gatherLocalIds = (iint*) calloc(nNonHalo, sizeof(iint));
  solver->nonHalo->gatherBaseIds  = (iint*) calloc(solver->nonHalo->Ngather, sizeof(iint));

  // only finds bases
  nHalo = 0;
  nNonHalo = 0;
  solver->halo->Ngather = 0; // reset counter
  solver->nonHalo->Ngather = 0; // reset counter

  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (mesh->gatherBaseIds[n] != mesh->gatherBaseIds[n-1]);

    // increment unique base counter and record index into shuffled list of nodes
    if(mesh->gatherHaloFlags[n]==1){
      if(test){
        solver->halo->gatherOffsets[solver->halo->Ngather] = nHalo;
        solver->halo->gatherBaseIds[solver->halo->Ngather] = mesh->gatherBaseIds[n];
        ++(solver->halo->Ngather);
      }
      solver->halo->gatherLocalIds[nHalo] = mesh->gatherLocalIds[n];
      ++nHalo;
    }
  }

  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (mesh->gatherBaseIds[n] != mesh->gatherBaseIds[n-1]);

    if(mesh->gatherHaloFlags[n]!=1){
      if(test){
        solver->nonHalo->gatherOffsets[solver->nonHalo->Ngather] = nNonHalo;
        ++(solver->nonHalo->Ngather);
      }
      solver->nonHalo->gatherLocalIds[nNonHalo] = mesh->gatherLocalIds[n];
      ++nNonHalo;
    }
  }
  solver->halo->gatherOffsets[solver->halo->Ngather] = nHalo;
  solver->nonHalo->gatherOffsets[solver->nonHalo->Ngather] = nNonHalo;

  // if there are halo nodes to gather
  if(solver->halo->Ngather){

    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(solver->halo->Ngather*sizeof(dfloat), NULL);
    solver->halo->gatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer(); // (char*) calloc(solver->halo->Ngather*sizeof(dfloat), sizeof(char));
    //printf("host gatherTmp = %p, sizeof(dfloat) = %d, Ngather = %d\n",  o_gatherTmpPinned.getMappedPointer(), sizeof(dfloat), solver->halo->Ngather);

    //    solver->halo->gatherTmp = (char*) calloc(solver->halo->Ngather*sizeof(dfloat), sizeof(char));

    solver->halo->o_gatherTmp      = mesh->device.malloc(solver->halo->Ngather*sizeof(dfloat),           solver->halo->gatherTmp);
    solver->halo->o_gatherOffsets  = mesh->device.malloc((solver->halo->Ngather+1)*sizeof(iint), solver->halo->gatherOffsets);
    solver->halo->o_gatherLocalIds = mesh->device.malloc(nHalo*sizeof(iint),                solver->halo->gatherLocalIds);

    // initiate gslib gather-scatter comm pattern on halo nodes only
    int verbose = strstr(options,"VERBOSE") ? 1:0;
    solver->halo->haloGsh = gsParallelGatherScatterSetup(solver->halo->Ngather, solver->halo->gatherBaseIds,verbose);
  }

  // if there are non-halo nodes to gather
  if(solver->nonHalo->Ngather){

    solver->nonHalo->haloGsh = NULL;

    solver->nonHalo->o_gatherOffsets  = mesh->device.malloc((solver->nonHalo->Ngather+1)*sizeof(iint), solver->nonHalo->gatherOffsets);
    solver->nonHalo->o_gatherLocalIds = mesh->device.malloc(nNonHalo*sizeof(iint),                solver->nonHalo->gatherLocalIds);
  }
  

  // count elements that contribute to global C0 gather-scatter
  iint globalCount = 0;
  iint localCount = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    int isHalo = 0;
    for(iint n=0;n<mesh->Np;++n){
      if(mesh->globalHaloFlags[e*mesh->Np+n]>0){
        isHalo = 1;
        break;
      }
    }
    globalCount += isHalo;
    localCount += 1-isHalo;
  }

  iint *globalGatherElementList    = (iint*) calloc(globalCount, sizeof(iint));
  iint *localGatherElementList = (iint*) calloc(localCount, sizeof(iint));

  globalCount = 0;
  localCount = 0;

  for(iint e=0;e<mesh->Nelements;++e){
    int isHalo = 0;
    for(iint n=0;n<mesh->Np;++n){
      if(mesh->globalHaloFlags[e*mesh->Np+n]>0){
        isHalo = 1;
        break;
      }
    }
    if(isHalo){
      globalGatherElementList[globalCount++] = e;
    } else{
      localGatherElementList[localCount++] = e;
    }
  }
  //printf("local = %d, global = %d\n", localCount, globalCount);

  solver->NglobalGatherElements = globalCount;
  solver->NlocalGatherElements = localCount;

  if(globalCount)
    solver->o_globalGatherElementList =
      mesh->device.malloc(globalCount*sizeof(iint), globalGatherElementList);

  if(localCount)
    solver->o_localGatherElementList =
      mesh->device.malloc(localCount*sizeof(iint), localGatherElementList);
}
