#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

void meshPrint3D(mesh3D *mesh){
  printf("EToV:\n");
  for(int e=0;e<mesh->Nelements;++e){
    printf("%d %d %d %d\n", 
	   mesh->EToV[e*mesh->Nverts+0],
	   mesh->EToV[e*mesh->Nverts+1],
	   mesh->EToV[e*mesh->Nverts+2],
	   mesh->EToV[e*mesh->Nverts+3]);
  }

  printf("EToE:\n");
  for(int e=0;e<mesh->Nelements;++e){
    printf("%d %d %d %d\n", 
	   mesh->EToE[e*mesh->Nfaces+0],
	   mesh->EToE[e*mesh->Nfaces+1],
	   mesh->EToE[e*mesh->Nfaces+2],
	   mesh->EToE[e*mesh->Nfaces+3]);
  }
}
