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

#include "elliptic.h"

int ellipticSolve(elliptic_t *elliptic, dfloat lambda, dfloat tol,
                  occa::memory &o_r, occa::memory &o_x){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int Niter = 0;
  int maxIter = 1000; 

  double start = 0.0, end =0.0;

#if 0
  if(options.compareArgs("VERBOSE","TRUE")){
    mesh->device.finish();
    start = MPI_Wtime(); 
  }
#endif
  
  Niter = pcg (elliptic, lambda, o_r, o_x, tol, maxIter);

#if 0
  if(options.compareArgs("VERBOSE","TRUE")){
    mesh->device.finish();
    end = MPI_Wtime();
    double localElapsed = end-start;

    if(mesh->rank==0) printf("Solver converged in %d iters \n", Niter );

    int size = mesh->size;

    hlong   localDofs = (hlong) mesh->Np*mesh->Nelements;
    hlong   localElements = (hlong) mesh->Nelements;
    double globalElapsed;
    hlong   globalDofs;
    hlong   globalElements;

    MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, mesh->comm );
    MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_HLONG,   MPI_SUM, 0, mesh->comm );
    MPI_Reduce(&localElements,&globalElements,1, MPI_HLONG,   MPI_SUM, 0, mesh->comm );

    if (mesh->rank==0){
      printf("%02d %02d "hlongFormat" "hlongFormat" %d %17.15lg %3.5g \t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME PRECONMEMORY] \n",
             mesh->size, mesh->N, globalElements, globalDofs, Niter, globalElapsed, elliptic->precon->preconBytes/(1E9));
    }
  }
#endif
  
  return Niter;

}
