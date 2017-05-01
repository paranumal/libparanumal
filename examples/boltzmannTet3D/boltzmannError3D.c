#include "boltzmann3D.h"

// currently maximum
void boltzmannError3D(mesh3D *mesh, dfloat time,char *options){
  dfloat maxQ1 = 0, minQ1 = 1e9;
  iint fid = 0; //  
  for(iint e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat q1=0;
      iint id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];

      maxQ1 = mymax(maxQ1, fabs(mesh->q[id*mesh->Nfields + fid]));
      minQ1 = mymin(minQ1, fabs(mesh->q[id*mesh->Nfields + fid]));
    }
  }

  // compute maximum over all processes
  dfloat globalMaxQ1, globalMinQ1;
  MPI_Allreduce(&maxQ1, &globalMaxQ1, 1,MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&minQ1, &globalMinQ1, 1,MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0)
    printf("%g %g %g (time,min(density),max(density)\n",time, globalMinQ1, globalMaxQ1);

  if(isnan(globalMinQ1) || isnan(globalMaxQ1))
    exit(EXIT_FAILURE);  
}
