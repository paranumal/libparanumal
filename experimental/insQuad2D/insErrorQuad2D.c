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

#include "insQuad2D.h"

// currently maximum
void insErrorQuad2D(ins_t *ins, dfloat time,char *options){

  mesh2D *mesh = ins->mesh;

  dfloat maxU = 0, minU = 1e9;
  dfloat maxV = 0, minV = 1e9;
  dfloat maxP = 0, minP = 1e9; 
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      dlong id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);
      maxU = mymax(maxU, fabs(ins->U[id]));
      minU = mymin(minU, fabs(ins->U[id]));
      //
      maxV = mymax(maxV, fabs(ins->V[id]));
      minV = mymin(minV, fabs(ins->V[id]));
      //
      maxP = mymax(maxP, fabs(ins->P[id]));
      minP = mymin(minP, fabs(ins->P[id]));
    }
  }

  // compute maximum over all processes
  dfloat gMaxU, gMinU;
  MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&minU, &gMinU, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  //
  dfloat gMaxV, gMinV;
  MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&minV, &gMinV, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  //
  dfloat gMaxP, gMinP;
  MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&minP, &gMinP, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0)
    printf("Time: %g minU: %g maxU: %g minV: %g maxV: %g minP: %g maxP: %g\n", 
             time, gMinU, gMaxU, gMinV, gMaxV,gMinP, gMaxP );

  if( isnan(gMinU) || isnan(gMaxU) || 
      isnan(gMinV) || isnan(gMaxV) || 
      isnan(gMinP) || isnan(gMaxP) ){
    occa::printTimer(); 
    exit(EXIT_FAILURE);
  }
}
