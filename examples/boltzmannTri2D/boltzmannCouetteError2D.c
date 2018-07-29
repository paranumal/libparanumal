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

#include "boltzmann2D.h"

// currently maximum
void boltzmannCouetteError2D(mesh2D *mesh, dfloat time){

  // Coutte Flow exact solution for U velocity

 dfloat maxerr = 0;
 dfloat maxQ1 = 0, minQ1 = 1e9;
 int fid = 1; //U velocity

dfloat nu = mesh->sqrtRT*mesh->sqrtRT/mesh->tauInv;

  for(int e=0;e<mesh->Nelements;++e)
  {
    for(int n=0;n<mesh->Np;++n)
    {
      dfloat q1=0;
      int id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      // U = sqrt(RT)*Q2/Q1; 
      dfloat u   = mesh->sqrtRT*mesh->q[id*mesh->Nfields + 1]/mesh->q[id*mesh->Nfields];

      
      dfloat uex = 2. * y ; 

        for(int k=1; k<=100; k++)
        {

         dfloat lamda = k*M_PI/0.5;
         uex += 2.*pow(-1,k)/(lamda*0.5)*exp(-nu*lamda*lamda*time)*sin(lamda*y);
        }

      maxerr = mymax(maxerr, fabs(u-uex));

      maxQ1 = mymax(maxQ1, fabs(mesh->q[id*mesh->Nfields]));
      minQ1 = mymin(minQ1, fabs(mesh->q[id*mesh->Nfields]));
      
    }
  }

  // compute maximum over all processes
  dfloat globalMaxQ1, globalMinQ1, globalMaxErr;
  MPI_Allreduce(&maxQ1, &globalMaxQ1, 1,
     MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&minQ1, &globalMinQ1, 1,
    MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxerr, &globalMaxErr, 1,
     MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);



  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    printf("%g %g %g %g (time,min(density),max(density),max(error)\n",
     time, globalMinQ1, globalMaxQ1, globalMaxErr);
    mesh->maxErrorBoltzmann = globalMaxErr; 
  }


  if(isnan(globalMaxErr))
    exit(EXIT_FAILURE);
  
}
