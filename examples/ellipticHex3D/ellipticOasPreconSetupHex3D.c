
#include "ellipticHex3D.h"

typedef struct{

  iint localId;
  iint baseId;
  iint maxRank;
  iint baseRank;

} preconGatherInfo_t;

int parallelCompareBaseRank(const void *a, const void *b){

  preconGatherInfo_t *fa = (preconGatherInfo_t*) a;
  preconGatherInfo_t *fb = (preconGatherInfo_t*) b;

  if(fa->baseRank < fb->baseRank) return -1;
  if(fa->baseRank > fb->baseRank) return +1;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;

}


ogs_t *ellipticOasPreconSetupHex3D(mesh3D *mesh){
  
  // assumes meshParallelGatherScatterSetup3D has been called
  
  // ????? need to extend storage for halo ?????
  
  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;
  
  // space for gather base indices with halo
  preconGatherInfo_t *gatherInfo = (preconGatherInfo_t*) calloc(Nlocal+Nhalo, sizeof(preconGatherInfo_t));

  // rearrange in node order
  for(iint n=0;n<Nlocal;++n){
    iint id = mesh->gatherLocalIds[n];
    gatherInfo[id].baseId   = mesh->gatherBaseIds[n];
    gatherInfo[id].maxRank  = mesh->gatherMaxRanks[n];
    gatherInfo[id].baseRank = mesh->gatherBaseRanks[n];
  }
  
  // exchange one element halo (fix later if warranted)
  preconGatherInfo_t *sendBuffer = (preconGatherInfo_t*) calloc(Nhalo, sizeof(preconGatherInfo_t));
  meshHaloExchange3D(mesh, mesh->Np*sizeof(preconGatherInfo_t), gatherInfo,   sendBuffer, gatherInfo+Nlocal);
  
  // offsetes to extract second layer
  iint *offset = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offset[0] = +mesh->Nfp;
  offset[1] = +mesh->Nq;
  offset[2] = -1;
  offset[3] = -mesh->Nq;
  offset[4] = +1;
  offset[5] = -mesh->Nfp;

  // build gather-scatter

  // a. gather on DEVICE   [ NuniqueBases on DEVICE, offsets, local node ids] 
  // b. get DEVICE=>DEVICE [ halo node ids ]
  // c. copy DEVICE=>HOST  
  // d. gs parallel gather scatter [ gsh ]
  // e. copy HOST=>DEVICE  
  // f. put on DEVICE      [ halo node ids ]
  // g. scatter on DEVICE  [ NuniqueBases on DEVICE, different than original gather ids ]
  
  // 1. create on-DEVICE gather arrays
  // 2. create on-DEVICE scatter arrays
  // 3. feed to meshParallelGatherScatterSetup3D.c (need to add scatter info Nscatter, ...)

  iint NqP = mesh->Nq+2;
  iint NpP = NqP*NqP*NqP;

  iint *offsetP = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offsetP[0] = +NqP*NqP;
  offsetP[1] = +NqP;
  offsetP[2] = -1;
  offsetP[3] = -NqP;
  offsetP[4] = +1;
  offsetP[5] = -NqP*NqP;

  iint *faceNodesPrecon = (iint*) calloc(mesh->Nfp*mesh->Nfaces, sizeof(iint));
  for(iint j=0;j<mesh->Nq;++j)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+j*mesh->Nq+0*mesh->Nfp] = i+1 + (j+1)*NqP + 0*NqP*NqP;
  
  for(iint k=0;k<mesh->Nq;++k)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+k*mesh->Nq+1*mesh->Nfp] = i+1 + 0*NqP + (k+1)*NqP*NqP;
  
  for(iint k=0;k<mesh->Nq;++k)
    for(iint j=0;j<mesh->Nq;++j)
      faceNodesPrecon[j+k*mesh->Nq+2*mesh->Nfp] = NqP-1 + (j+1)*NqP + (k+1)*NqP*NqP;

  for(iint k=0;k<mesh->Nq;++k)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+k*mesh->Nq+3*mesh->Nfp] = i+1 + (NqP-1)*NqP + (k+1)*NqP*NqP;

  for(iint k=0;k<mesh->Nq;++k)
    for(iint j=0;j<mesh->Nq;++j)
      faceNodesPrecon[j+k*mesh->Nq+4*mesh->Nfp] = 0 + (j+1)*NqP + (k+1)*NqP*NqP;
  
  for(iint j=0;j<mesh->Nq;++j)
    for(iint i=0;i<mesh->Nq;++i)
      faceNodesPrecon[i+j*mesh->Nq+5*mesh->Nfp] = i+1 + (j+1)*NqP + (NqP-1)*NqP*NqP;

  preconGatherInfo_t *preconGatherInfo = (preconGatherInfo_t*) calloc(NpP*mesh->Nelements, sizeof(preconGatherInfo_t));
  
  // 0 numbering for uninvolved nodes
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint k=0;k<mesh->Nq;++k){
      for(iint j=0;j<mesh->Nq;++j){
	for(iint i=0;i<mesh->Nq;++i){
	  iint id  = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq + e*mesh->Np;
	  iint pid = i + 1 + (j+1)*NqP + (k+1)*NqP*NqP + e*NpP;

	  // ugly - need to make a struct then sort
	  preconGatherInfo[pid] = gatherInfo[id];
	}
      }
    }
  }

  // take node info from positive trace and put in overlap region on parallel gather info
  for(iint e=0;e<mesh->Nelements;++e){

    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      iint fid = e*mesh->Nfp*mesh->Nfaces + n;
      iint fM = n/mesh->Nfp;
      iint fP = mesh->EToF[e*mesh->Nfaces+fM];
      if(fP<0) fP = fM;

      iint idP = mesh->vmapP[fid] + offset[fP];
      iint idO = faceNodesPrecon[n] + e*NpP;
      preconGatherInfo[idO] = gatherInfo[idP];
    }
  }
  
#if 1
  for(iint e=0;e<mesh->Nelements;++e){
    printf("e=%d:[ \n", e);
    for(iint k=0;k<NqP;++k){
      printf("     [\n");
      for(iint j=0;j<NqP;++j){
	printf("       ");
	for(iint i=0;i<NqP;++i){
	  iint id  = i + j*NqP + k*NqP*NqP + e*NpP;
	  printf("%05d ", preconGatherInfo[id].baseId);
	}
	printf("\n");
      }
      printf("     ]\n");
    }
    printf(" ]\n");
  }
#endif
  
  // reset local ids
  for(iint n=0;n<mesh->Nelements*NpP;++n)
    preconGatherInfo[n].localId = n;
  
  // sort by rank then base index
  qsort(preconGatherInfo, NpP*mesh->Nelements, sizeof(preconGatherInfo_t), parallelCompareBaseRank);

  // do not gather-scatter nodes labelled zero
  iint skip = 0;
  while(preconGatherInfo[skip].baseId==0 && skip<NpP*mesh->Nelements){
    ++skip;
  }

  printf("skip = %d, NlocalP = %d\n", skip, NpP*mesh->Nelements);
  
  // reset local ids
  iint NlocalP = NpP*mesh->Nelements - skip;
  iint *gatherLocalIdsP  = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherBaseIdsP   = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherBaseRanksP = (iint*) calloc(NlocalP, sizeof(iint));
  iint *gatherMaxRanksP  = (iint*) calloc(NlocalP, sizeof(iint));
  for(iint n=0;n<NlocalP;++n){
    gatherLocalIdsP[n]  = preconGatherInfo[n+skip].localId;
    gatherBaseIdsP[n]   = preconGatherInfo[n+skip].baseId;
    gatherBaseRanksP[n] = preconGatherInfo[n+skip].baseRank;
    gatherMaxRanksP[n]  = preconGatherInfo[n+skip].maxRank;
    //    printf("base[%d] = %d\n", n, gatherBaseIdsP[n]);
  }
  
  // make preconBaseIds => preconNumbering
  ogs_t *ogsP = meshParallelGatherScatterSetup3D(mesh, NlocalP, sizeof(dfloat),
						 gatherLocalIdsP, gatherBaseIdsP, gatherBaseRanksP, gatherMaxRanksP);

  // as is - will need to extract local nodes from precon nodes
  
  return ogsP;
}
