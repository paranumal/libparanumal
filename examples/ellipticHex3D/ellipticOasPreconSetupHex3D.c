
#include "ellipticHex3D.h"

typedef struct{

  iint localId;
  iint baseId;
  iint maxRank;
  iint baseRank;

} gatherInfo_t;

int parallelCompareBaseRank(const void *a, const void *b){

  gatherInfo_t *fa = (gatherInfo_t*) a;
  gatherInfo_t *fb = (gatherInfo_t*) b;

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
  gatherInfo_t *gatherInfo = (gatherInfo_t*) calloc(Nlocal+Nhalo, sizeof(gatherInfo_t));

  // rearrange in node order
  for(iint n=0;n<Nlocal;++n){
    iint id = mesh->gatherLocalIds[n];
    gatherInfo[id].baseId   = mesh->gatherBaseIds[n];
    gatherInfo[id].maxRank  = mesh->gatherMaxRanks[n];
    gatherInfo[id].baseRank = mesh->gatherBaseRanks[n];
  }
  
  // exchange one element halo (fix later if warranted)
  gatherInfo_t *sendBuffer = (gatherInfo_t*) calloc(Nhalo, sizeof(gatherInfo_t));
  meshHaloExchange3D(mesh, mesh->Np*sizeof(gatherInfo_t), gatherInfo,   sendBuffer, gatherInfo+Nlocal);
  
  // offsetes to extract second layer
  iint *offset = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offset[0] = +mesh->Nfp;
  offset[1] = +mesh->Nq;
  offset[2] = -1;
  offset[3] = -mesh->Nq;
  offset[4] = +1;
  offset[5] = -mesh->Nfp;

#if 0
  iint *offsetP = (iint*) calloc(mesh->Nfaces, sizeof(iint));
  offsetP[0] = +NqP*NqP;
  offsetP[1] = +NqP;
  offsetP[2] = -1;
  offsetP[3] = -NqP;
  offsetP[4] = +1;
  offsetP[5] = -NqP*NqP;
#endif
  
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

#if 0
  // need to build these for gather and scatter steps
  iint Ngather;
  iint *gatherLocalIds;  // local index of nodes
  iint *gatherBaseIds;   // global index of their base nodes
  iint *gatherBaseRanks; // rank of their base nodes
  iint *gatherMaxRanks;  // max rank connected to base node

  // HMMM CAN JUST USE THESE AND THEN EXTRACT THE INTERIOR CHUNKS
#endif
  
  iint NqP = mesh->Nq+2;
  iint NpP = NqP*NqP*NqP;
  gatherInfo_t *gatherInfoP = (gatherInfo_t*) calloc(NpP*mesh->Nelements, sizeof(gatherInfo_t));

  // 0 numbering for uninvolved nodes
  for(iint e=0;e<mesh->Nelements;++e){
    
    for(iint k=0;k<mesh->Nq;++k){
      for(iint j=0;j<mesh->Nq;++j){
	for(iint i=0;i<mesh->Nq;++i){
	  iint id  = i + j*mesh->Nq + k*mesh->Nq*mesh->Nq + e*mesh->Np;
	  iint pid = i + 1 + (j+1)*NqP + (k+1)*NqP*NqP + e*NpP;

	  // ugly - need to make a struct then sort
	  gatherInfoP[pid] = gatherInfo[id];

	  // brittle
	  if(k==0){
	    iint fid = i + j*mesh->Nq + 0*mesh->Nfp + e*mesh->Nfp*mesh->Nfaces;

	    iint fP = mesh->EToF[e*mesh->Nfaces+0];
	    if(fP<0) fP = 0;

	    iint idP =  mesh->vmapP[fid] + offset[fP]; // need to use offset of plus face
	    gatherInfoP[pid-NqP*NqP] = gatherInfo[idP];  // different offsets
	  }
	  if(j==0){
	    iint fid = i + k*mesh->Nq + 1*mesh->Nfp + e*mesh->Nfp*mesh->Nfaces;

	    iint fP = mesh->EToF[e*mesh->Nfaces+1];
	    if(fP<0) fP = 1;

	    iint idP = mesh->vmapP[fid] + offset[fP];
	    gatherInfoP[pid - NqP] = gatherInfo[idP];
	  }
	  if(i==mesh->Nq-1){
	    iint fid = j + k*mesh->Nq + 2*mesh->Nfp + e*mesh->Nfp*mesh->Nfaces;

	    iint fP = mesh->EToF[e*mesh->Nfaces+2];
	    if(fP<0) fP = 2;

	    iint idP = mesh->vmapP[fid] + offset[fP];
	    gatherInfoP[pid + 1] = gatherInfo[idP];
	  }
	  if(j==mesh->Nq-1){
	    iint fid = i + k*mesh->Nq + 3*mesh->Nfp + e*mesh->Nfp*mesh->Nfaces;

	    iint fP = mesh->EToF[e*mesh->Nfaces+3];
	    if(fP<0) fP = 3;

	    iint idP = mesh->vmapP[fid] + offset[fP];
	    gatherInfoP[pid + NqP] = gatherInfo[idP];
	  }
	  if(i==0){
	    iint fid = j + k*mesh->Nq + 4*mesh->Nfp + e*mesh->Nfp*mesh->Nfaces;

	    iint fP = mesh->EToF[e*mesh->Nfaces+4];
	    if(fP<0) fP = 4;

	    iint idP = mesh->vmapP[fid] + offset[fP];
	    gatherInfoP[pid - 1] = gatherInfo[idP];
	  }
	  if(k==mesh->Nq-1){
	    iint fid = i + j*mesh->Nq + 5*mesh->Nfp + e*mesh->Nfp*mesh->Nfaces;	    

	    iint fP = mesh->EToF[e*mesh->Nfaces+5];
	    if(fP<0) fP = 5;

	    iint idP = mesh->vmapP[fid] + offset[fP];
	    //	    printf("mesh->vmapP[fid] = %d\n", mesh->vmapP[fid]);
	    gatherInfoP[pid + NqP*NqP] = gatherInfo[idP];
	  }
	}
      }
    }
  }
#if 0
  for(iint e=0;e<mesh->Nelements;++e){
    printf("e=%d:[ \n", e);
    for(iint k=0;k<NqP;++k){
      printf("     [\n");
      for(iint j=0;j<NqP;++j){
	printf("       ");
	for(iint i=0;i<NqP;++i){
	  iint id  = i + j*NqP + k*NqP*NqP + e*NpP;
	  printf("%05d ", gatherInfoP[id].baseId);
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
    gatherInfoP[n].localId = n;
  
  // sort by rank then base index
  qsort(gatherInfoP, NpP*mesh->Nelements, sizeof(gatherInfo_t), parallelCompareBaseRank);

  // do not gather-scatter nodes labelled zero
  iint skip = 0;
  while(gatherInfoP[skip].baseId==0 && skip<NpP*mesh->Nelements){
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
    gatherLocalIdsP[n]  = gatherInfoP[n+skip].localId;
    gatherBaseIdsP[n]   = gatherInfoP[n+skip].baseId;
    gatherBaseRanksP[n] = gatherInfoP[n+skip].baseRank;
    gatherMaxRanksP[n]  = gatherInfoP[n+skip].maxRank;
    //    printf("base[%d] = %d\n", n, gatherBaseIdsP[n]);
  }
  
  // make preconBaseIds => preconNumbering
  ogs_t *ogsP = meshParallelGatherScatterSetup3D(mesh, NlocalP, sizeof(dfloat),
						 gatherLocalIdsP, gatherBaseIdsP, gatherBaseRanksP, gatherMaxRanksP);

  // as is - will need to extract local nodes from precon nodes
  
  return ogsP;
}
