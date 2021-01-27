/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondAMGSetup.hpp"

namespace parAlmond {

void parAlmond_t::AMGSetup(parCOO& cooA,
                         bool nullSpace,
                         dfloat *nullVector,
                         dfloat nullSpacePenalty){

  int rank;
  int size;
  MPI_Comm_rank(cooA.comm, &rank);
  MPI_Comm_size(cooA.comm, &size);

  if(rank==0) {printf("Setting up AMG...");fflush(stdout);}

  //make csr matrix from coo input
  parCSR *A = new parCSR(cooA);
  A->diagSetup();

  //copy fine nullvector
  dfloat *null = (dfloat *) malloc(A->Nrows*sizeof(dfloat));
  memcpy(null, nullVector, A->Nrows*sizeof(dfloat));

  // find target N at coarsest level
  const int gCoarseSize = multigrid->coarseSolver->getTargetSize();

  amgLevel *L = new amgLevel(A, settings);

  hlong globalSize;
  if (multigrid->coarsetype==COARSEEXACT) {
    globalSize = L->A->globalRowStarts[size];
  } else { //COARSEOAS
    //OAS cares about Ncols for size
    hlong localSize = A->Ncols;
    MPI_Allreduce(&localSize,&globalSize,1,MPI_HLONG,MPI_SUM,A->comm);
  }

  //if the system if already small, dont create MG levels
  bool done = false;
  if(globalSize <= gCoarseSize){
    multigrid->AddLevel(L);
    multigrid->coarseSolver->setup(A, nullSpace, null, nullSpacePenalty);
    multigrid->coarseSolver->syncToDevice();
    multigrid->baseLevel = multigrid->numLevels-1;
    L->syncToDevice();
    done = true;
  }

  //TODO: make the coarsen threasholds user-provided inputs
  // For now, let default to some sensible threasholds
  dfloat theta=0.0;
  if (multigrid->strtype==RUGESTUBEN) {
    //    theta=0.14; //default for 3D problems
    //    theta=0.12; //default for 3D problems
    theta = 0.5;
    settings.getSetting("PARALMOND RUGESTUBEN STRENGTH THRESHOLD", theta);
    printf("theta = %f\n", theta);
    //See: A GPU accelerated aggregation algebraic multigrid method, R. Gandham, K. Esler, Y. Zhang.
  } else { // (type==SYMMETRIC)
    theta=0.08;
    //See: Algebraic Multigrid On Unstructured Meshes, P Vanek, J. Mandel, M. Brezina.
  }

  while(!done){
    L->setupSmoother();

    // Create coarse level via AMG. Coarsen null vector
    amgLevel* Lcoarse = coarsenAmgLevel(L, null,
                                        multigrid->strtype, theta,
                                        multigrid->aggtype);
    multigrid->AddLevel(L);
    L->syncToDevice();

    // Increase coarsening rate as we add levels.
    //See: Algebraic Multigrid On Unstructured Meshes, P Vanek, J. Mandel, M. Brezina.
    if (multigrid->strtype==SYMMETRIC)
      theta=theta/2;

    hlong globalCoarseSize;
    if (multigrid->coarsetype==COARSEEXACT) {
      globalCoarseSize = Lcoarse->A->globalRowStarts[size];;
    } else { //COARSEOAS
      //OAS cares about Ncols for size
      hlong localSize = Lcoarse->A->Ncols;
      MPI_Allreduce(&localSize,&globalCoarseSize,1,MPI_HLONG,MPI_SUM,Lcoarse->A->comm);
    }

    if(globalCoarseSize <= gCoarseSize || globalSize < 2*globalCoarseSize){
      if (globalSize < 2*globalCoarseSize && rank==0) {
        stringstream ss;
        ss << "AMG coarsening stalling, attemping coarse solver setup with dimension N=" << globalCoarseSize;
        LIBP_WARNING(ss.str());
      }
      multigrid->AddLevel(Lcoarse);
      Lcoarse->syncToDevice();
      multigrid->coarseSolver->setup(Lcoarse->A, nullSpace, null, nullSpacePenalty);
      multigrid->coarseSolver->syncToDevice();
      multigrid->baseLevel = multigrid->numLevels-1;
      break;
    }
    globalSize = globalCoarseSize;
    L = Lcoarse;
  }

  free(null);

  if(rank==0) printf("done.\n");
}

} //namespace parAlmond
