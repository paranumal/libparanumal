#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

void meshPrint3D(mesh3D *mesh){
  printf("EToV:\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int v=0;v<mesh->Nverts;++v){
      printf(hlongFormat " ", mesh->EToV[e*mesh->Nverts+v]);
    }
    printf("\n");
  }

  printf("EToE:\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      printf(dlongFormat " ",  mesh->EToE[e*mesh->Nfaces+f]);
    }
    printf("\n");
  }

  printf("EToB:\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      printf("%d ",  mesh->EToB[e*mesh->Nfaces+f]);
    }
    printf("\n");
  }
}
