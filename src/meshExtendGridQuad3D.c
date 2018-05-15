#include "mesh.h"

void meshExtendGridQuad3D(mesh3D *mesh) {
  mesh->NgridElements = mesh->Nelements + mesh->edgeLength*24;
  mesh->gridToE = (iint *) calloc(mesh->Nelements*mesh->Nfaces,sizeof(iint));
  mesh->overlap = (iint *) calloc((mesh->NgridElements - mesh->Nelements),sizeof(iint));
  iint count_extra = 0;
  for (iint e = 0; e < mesh->Nelements;++e) {
    for (iint f = 0; f < mesh->Nfaces; ++f) {
      if (mesh->cubeFaceNumber[e] == mesh->cubeFaceNumber[mesh->EToE[e*mesh->Nfaces + f]]) {
	mesh->gridToE[e*mesh->Nfaces + f] = mesh->EToE[e*mesh->Nfaces + f];
      }
      else {
	mesh->gridToE[e*mesh->Nfaces + f] = mesh->Nelements + count_extra;
	mesh->overlap[count_extra] = mesh->cubeFaceNumber[mesh->EToE[e*mesh->Nfaces + f]];
	count_extra++;
      }
  }
}

void meshPreserveGridQuad3D(mesh3D *mesh) {
  mesh->NgridElements = mesh->Nelements + mesh->edgeLength*24;
  mesh->gridToE = mesh->EToE; //shared pointer
  //count_extra is not needed
}
