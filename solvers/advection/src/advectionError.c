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

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "advection.h"


void advectionError(advection_t *advection, dfloat time){

  mesh_t *mesh = advection->mesh;

  dfloat maxR = 0;
  dfloat minR = 1E9;
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dfloat u,v,w, p;
      int id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];

      int qbase = n+e*mesh->Np*mesh->Nfields;
      maxR = mymax(maxR, advection->q[qbase]);
      minR = mymin(minR, advection->q[qbase]);
    }
  }

  // compute maximum over all processes
  dfloat globalMaxR;
  dfloat globalMinR;
  MPI_Allreduce(&maxR, &globalMaxR, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
  MPI_Allreduce(&minR, &globalMinR, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

  if(mesh->rank==0)
    printf("%g, %g, %g ( time, min density, max density)\n", time, globalMinR, globalMaxR);
  
}
