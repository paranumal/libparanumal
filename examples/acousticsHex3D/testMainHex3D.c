#include "acousticsHex3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // set up mesh stuff
  mesh3D *meshSetupHex3D(char *, int);
  mesh3D *mesh = meshSetupHex3D(argv[1], N);

  // set up acoustics stuff
  void acousticsSetupHex3D(mesh3D *mesh);
  acousticsSetupHex3D(mesh);

  dfloat *p = (dfloat*) calloc(mesh->Np*mesh->Nelements, sizeof(dfloat));

  dfloat sumError = 0;
  
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n)
      p[e*mesh->Np+n] = mesh->z[e*mesh->Np+n];

    for(int n=0;n<mesh->Np;++n){
      dfloat dpdr = 0, dpds = 0, dpdt = 0;
      for(int m=0;m<mesh->Np;++m){
	dpdr += mesh->Dr[n*mesh->Np+m]*p[e*mesh->Np+m];
	dpds += mesh->Ds[n*mesh->Np+m]*p[e*mesh->Np+m];
	dpdt += mesh->Dt[n*mesh->Np+m]*p[e*mesh->Np+m];
      }
      int base = e*mesh->Np*mesh->Nvgeo + n;
      dfloat drdx = mesh->vgeo[base + RXID*mesh->Np];
      dfloat dsdx = mesh->vgeo[base + SXID*mesh->Np];
      dfloat dtdx = mesh->vgeo[base + TXID*mesh->Np];
      dfloat drdy = mesh->vgeo[base + RYID*mesh->Np];
      dfloat dsdy = mesh->vgeo[base + SYID*mesh->Np];
      dfloat dtdy = mesh->vgeo[base + TYID*mesh->Np];
      dfloat drdz = mesh->vgeo[base + RZID*mesh->Np];
      dfloat dsdz = mesh->vgeo[base + SZID*mesh->Np];
      dfloat dtdz = mesh->vgeo[base + TZID*mesh->Np];
      dfloat dpdx = drdx*dpdr + dsdx*dpds + dtdx*dpdt;
      dfloat dpdy = drdy*dpdr + dsdy*dpds + dtdy*dpdt;
      dfloat dpdz = drdz*dpdr + dsdz*dpds + dtdz*dpdt;
      //      printf("dpdx = %g, dpdy = %g, dpdz = %g\n",
      //	     dpdx, dpdy, dpdz);
      sumError += (dpdz-1)*(dpdz-1);
      //      if(isnan(dpdz)) printf("nan-nan-nan\n");
    }
  }

  printf("L2 error = %g\n", sqrt(sumError));
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
