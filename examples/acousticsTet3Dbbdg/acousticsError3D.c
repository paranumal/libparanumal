#include "acoustics3D.h"

void acousticsError3D(mesh3D *mesh, dfloat time){

  dfloat maxErrorP = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      int id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];

      dfloat u,v,w,p;
      //acousticsCavitySolution3D(x, y, z, time, &u, &v, &w, &p);
      p = 0.f;

      dfloat qm = mesh->q[(e*mesh->Np+n)*mesh->Nfields+3];

      maxErrorP = mymax(maxErrorP, fabs(p-qm));
    }
  }

  // compute maximum over all processes
  dfloat globalMaxErrorP;
  MPI_Allreduce(&maxErrorP, &globalMaxErrorP, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0)
    printf("%g, %g (time,maxError(pressure)\n", time, globalMaxErrorP);
  
}
