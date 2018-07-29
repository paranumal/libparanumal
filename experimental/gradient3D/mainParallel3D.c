/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

void message(char *s){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("%s", s);
  MPI_Barrier(MPI_COMM_WORLD);
}


int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  if(argc!=3){
    printf("usage: ./main meshname.msh N\n");
    exit(-1);
  }

  mesh3D *mesh;

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  mesh = meshSetup3D(argv[1], N);

  // compute samples of q at interpolation nodes
  dfloat *q = (dfloat*) calloc(mesh->Nelements*mesh->Np,
			       sizeof(dfloat));
  dfloat *dqdx = (dfloat*) calloc(mesh->Nelements*mesh->Np,
				  sizeof(dfloat));
  dfloat *dqdy = (dfloat*) calloc(mesh->Nelements*mesh->Np,
				  sizeof(dfloat));
  dfloat *dqdz = (dfloat*) calloc(mesh->Nelements*mesh->Np,
				  sizeof(dfloat));
  
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      q[cnt] = mesh->z[cnt];
      ++cnt;
    }
  }

  // compute interpolant and collocation derivative of
  // physical x,y derivatives of q data 
  // interpolated at interpolation nodes
  meshGradient3D(mesh, q, dqdx, dqdy, dqdz);

  // run a test kernel
  occaTest3D(mesh, q, dqdx, dqdy, dqdz);

  // run timings for 9 versions of meshGradient kernel 
  occaOptimizeGradient3D(mesh, q, dqdx, dqdy, dqdz);

  // estimate pointwise error at interpolation nodes
  dfloat maxError = 0;
  cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      maxError = mymax(maxError, fabs(dqdz[cnt]-1));
      ++cnt;
    }
  }

  printf("maxError = %6.5E\n", maxError);

  // create single vtu file name
  char vtuName[BUFSIZ];
  sprintf(vtuName, "%s.vtu", strtok(argv[1], "."));
  meshVTU3D(mesh, vtuName);
  
  MPI_Finalize();

  exit(0);
  return 0;
}
