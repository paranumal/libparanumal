
#include "ellipticHex3D.h"

void ellipticOasPreconSetupHex3D(mesh3D *mesh){

  // ????? need to extend storage for halo ?????
  
  iint Nlocal = mesh->Np*mesh->Nelements;
  iint Nhalo  = mesh->Np*mesh->totalHaloPairs;
  iint Ntrace = mesh->Nfp*mesh->Nfaces*mesh->Nelements;
  
  // space for global base indices with halo
  iint *globalBaseIds = (iint*) calloc(Nlocal+Nhalo, sizeof(iint));

  // populate local numbers
  memcpy(globalBaseIds, mesh->globalBaseIds, Nlocal*sizeof(iint));

  // exchange DG halo
  iint *sendBuffer = (iint*) calloc(Nhalo, sizeof(iint));
  meshHaloExchange3D(mesh, mesh->Np*sizeof(iint), globalBaseIds, sendBuffer,
		     globalBaseIds+Nlocal);
  
  // extract second layer
  iint *offset2 = (iint*) calloc(mesh->Nfaces*mesh->Nfp, sizeof(iint));
  for(iint n=0;n<mesh->Nfp;++n) offset2[n + 0*mesh->Nfp] = +mesh->Nfp;
  for(iint n=0;n<mesh->Nfp;++n) offset2[n + 1*mesh->Nfp] = +mesh->Nq;
  for(iint n=0;n<mesh->Nfp;++n) offset2[n + 2*mesh->Nfp] = -1;
  for(iint n=0;n<mesh->Nfp;++n) offset2[n + 3*mesh->Nfp] = -mesh->Nq;
  for(iint n=0;n<mesh->Nfp;++n) offset2[n + 4*mesh->Nfp] = +1;
  for(iint n=0;n<mesh->Nfp;++n) offset2[n + 5*mesh->Nfp] = -mesh->Nfp;
  
  // collect vmapPP from globalBaseIds
  mesh->vmapPP = (iint*) calloc(Ntrace, sizeof(iint));

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      iint id = e*mesh->Nfp*mesh->Nfaces + n;
      iint idP = mesh->vmapP[id];
      mesh->vmapPP[id] = globalBaseIds[idP+offset2[n]];
    }
  }

  mesh->o_vmapPP = mesh->device.malloc(Ntrace*sizeof(iint), mesh->vmapPP);

  //      x a a a a x
  //      b e e e e c
  //      b e e e e c
  //      b e e e e c
  //      b e e e e c
  //      x d d d d x
  
  // build gather-scatter

  // a. gather on DEVICE   [ NuniqueBases on DEVICE, offsets, local node ids] 
  // b. get DEVICE=>DEVICE [ halo node ids ]
  // c. copy DEVICE=>HOST  
  // d. gs parallel gather scatter [ gsh ]
  // e. copy HOST=>DEVICE  
  // f. put on DEVICE      [ halo node ids ]
  // g. scatter on DEVICE  [ NuniqueBases on DEVICE, different than original gather ids ]
  
  // 1. create on-DEVICE gather arrays
  // 2. 
  
}
