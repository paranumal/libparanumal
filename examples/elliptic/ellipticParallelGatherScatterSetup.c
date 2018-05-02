#include "elliptic.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(elliptic_t* elliptic){  

  mesh_t *mesh = elliptic->mesh;

  // setup occa gather scatter
  int verbose = elliptic->options.compareAgrs("VERBOSE","TRUE") ? 1:0;
  mesh->ogs = meshParallelGatherScatterSetup(mesh,mesh->Np*mesh->Nelements,
                                             mesh->gatherLocalIds,
                                             mesh->gatherBaseIds,
                                             mesh->gatherBaseRanks,
                                             mesh->gatherHaloFlags,
                                             verbose);
  elliptic->o_invDegree = mesh->ogs->o_invDegree;

  // count elements that contribute to global C0 gather-scatter
  dlong globalCount = 0;
  dlong localCount = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
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

  dlong *globalGatherElementList = (dlong*) calloc(globalCount, sizeof(dlong));
  dlong *localGatherElementList  = (dlong*) calloc(localCount, sizeof(dlong));

  globalCount = 0;
  localCount = 0;

  for(dlong e=0;e<mesh->Nelements;++e){
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

  elliptic->NglobalGatherElements = globalCount;
  elliptic->NlocalGatherElements = localCount;

  if(globalCount)
    elliptic->o_globalGatherElementList =
      mesh->device.malloc(globalCount*sizeof(dlong), globalGatherElementList);

  if(localCount)
    elliptic->o_localGatherElementList =
      mesh->device.malloc(localCount*sizeof(dlong), localGatherElementList);
}
