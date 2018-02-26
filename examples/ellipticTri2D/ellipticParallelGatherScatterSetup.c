#include "ellipticTri2D.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(solver_t* solver, const char *options){  

  mesh2D *mesh = solver->mesh;

  // setup occa gather scatter
  int verbose = strstr(options,"VERBOSE") ? 1:0;
  mesh->ogs = meshParallelGatherScatterSetup(mesh,mesh->Np*mesh->Nelements,
                                             mesh->gatherLocalIds,
                                             mesh->gatherBaseIds,
                                             mesh->gatherBaseRanks,
                                             mesh->gatherHaloFlags,
                                             verbose);
  solver->o_invDegree = mesh->ogs->o_invDegree;

  // count elements that contribute to global C0 gather-scatter
  int globalCount = 0;
  int localCount = 0;
  for(int e=0;e<mesh->Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<mesh->Np;++n){
      if(mesh->globalHaloFlags[e*mesh->Np+n]>0){
        isHalo = 1;
        break;
      }
    }
    globalCount += isHalo;
    localCount += 1-isHalo;
  }

  int *globalGatherElementList    = (int*) calloc(globalCount, sizeof(int));
  int *localGatherElementList = (int*) calloc(localCount, sizeof(int));

  globalCount = 0;
  localCount = 0;

  for(int e=0;e<mesh->Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<mesh->Np;++n){
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
      mesh->device.malloc(globalCount*sizeof(int), globalGatherElementList);

  if(localCount)
    solver->o_localGatherElementList =
      mesh->device.malloc(localCount*sizeof(int), localGatherElementList);
}
