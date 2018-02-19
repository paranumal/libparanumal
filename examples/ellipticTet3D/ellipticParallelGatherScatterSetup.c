#include "ellipticTet3D.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(solver_t* solver, const char *options){  

  mesh3D *mesh = solver->mesh;

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
