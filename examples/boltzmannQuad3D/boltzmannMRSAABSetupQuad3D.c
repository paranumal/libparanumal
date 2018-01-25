#include "mesh3D.h"

//break out setup needed by the multirate code
void boltzmannMRSAABSetupQuad3D (mesh_t *mesh) {
  dfloat *EtoDT       = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));

  for (iint e = 0; e < mesh->Nelements; ++e) {
    EtoDT[e] = 1e9; //dummy initial value
    dfloat hmin = 1e9;

    for(iint f=0;f<mesh->Nfaces;++f){
      iint sid    = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
      
      dfloat htest = 0.5/(sJ*invJ);
      hmin = mymin(hmin, htest); 
    }
  
  meshMRABSetup3D(mesh, 
  
  mesh->MRAB_A    = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  mesh->MRSAAB_A  = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  mesh->MRSAAB_C  = (dfloat *) calloc(    Nlevels,sizeof(dfloat));
