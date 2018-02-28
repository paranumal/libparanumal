#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

void meshPrint3D(mesh3D *mesh){
  printf("EToV:\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    printf(hlongFormat hlongFormat hlongFormat hlongFormat"\n", 
	   mesh->EToV[e*mesh->Nverts+0],
	   mesh->EToV[e*mesh->Nverts+1],
	   mesh->EToV[e*mesh->Nverts+2],
	   mesh->EToV[e*mesh->Nverts+3]);
  }

  printf("EToE:\n");
  for(dlong e=0;e<mesh->Nelements;++e){
    printf(dlongFormat dlongFormat dlongFormat dlongFormat"\n", 
	   mesh->EToE[e*mesh->Nfaces+0],
	   mesh->EToE[e*mesh->Nfaces+1],
	   mesh->EToE[e*mesh->Nfaces+2],
	   mesh->EToE[e*mesh->Nfaces+3]);
  }
}
