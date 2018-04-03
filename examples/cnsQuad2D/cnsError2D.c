#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"


void cnsError2D(mesh2D *mesh, dfloat time){

  dfloat maxErrorP = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat u,v,p;
      int id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      int qbase = n+e*mesh->Np*mesh->Nfields;
      u = 0;
      v = 0;
      p = 0;

      maxErrorP = mymax(maxErrorP, fabs(p-mesh->q[qbase]));
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
