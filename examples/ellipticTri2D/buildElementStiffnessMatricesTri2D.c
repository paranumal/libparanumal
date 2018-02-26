#include "ellipticTri2D.h"

void buildElementStiffnessMatricesTri2D(mesh2D *mesh, const char *options, int N){



  dfloat *SrrT, *SrsT, *SsrT, *SssT;
  mesh->Srr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  mesh->Srs = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  mesh->Ssr = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  mesh->Sss = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (iint n=0;n<mesh->Np;n++) {
    for (iint m=0;m<mesh->Np;m++) {
      for (iint k=0;k<mesh->Np;k++) {
        for (iint l=0;l<mesh->Np;l++) {
          mesh->Srr[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
          mesh->Srs[m+n*mesh->Np] += mesh->Dr[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
          mesh->Ssr[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Dr[m+k*mesh->Np];
          mesh->Sss[m+n*mesh->Np] += mesh->Ds[n+l*mesh->Np]*mesh->MM[k+l*mesh->Np]*mesh->Ds[m+k*mesh->Np];
        }
      } 
    }
  }
  SrrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  SrsT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  SsrT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  SssT = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));
  for (iint n=0;n<mesh->Np;n++) {
    for (iint m=0;m<mesh->Np;m++) {  
      SrrT[m+n*mesh->Np] = mesh->Srr[n+m*mesh->Np];
      SrsT[m+n*mesh->Np] = mesh->Srs[n+m*mesh->Np];
      SsrT[m+n*mesh->Np] = mesh->Ssr[n+m*mesh->Np];
      SssT[m+n*mesh->Np] = mesh->Sss[n+m*mesh->Np];
      //KS
    }
  }
  //

  mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrrT);
  mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SrsT);
  mesh->o_SsrT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SsrT);
  mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->Np*sizeof(dfloat), SssT);

}



