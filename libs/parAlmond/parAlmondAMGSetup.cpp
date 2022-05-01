/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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
#include "parAlmond/parAlmondCoarseSolver.hpp"

namespace libp {

namespace parAlmond {

void parAlmond_t::AMGSetup(parCOO& cooA,
                         bool nullSpace,
                         memory<dfloat> nullVector,
                         dfloat nullSpacePenalty){

  int rank = cooA.comm.rank();
  int size = cooA.comm.size();

  if(Comm::World().rank()==0) {printf("Setting up AMG...");fflush(stdout);}

  /*Get multigrid solver*/
  multigrid_t& mg = *multigrid;

  /*Get coarse solver*/
  coarseSolver_t& coarse = *(mg.coarseSolver);

  //make csr matrix from coo input
  parCSR A(cooA);
  A.diagSetup();

  //copy fine nullvector
  memory<dfloat> null(A.Nrows);
  null.copyFrom(nullVector, A.Nrows);

  // find target N at coarsest level
  const int gCoarseSize = coarse.getTargetSize();

  hlong globalSize;
  if (mg.coarsetype==COARSEEXACT) {
    globalSize = A.globalRowStarts[size];
  } else { //COARSEOAS
    //OAS cares about Ncols for size
    globalSize = A.Ncols;
    A.comm.Allreduce(globalSize);
  }

  amgLevel& Lbase = mg.AddLevel<amgLevel>(A, settings);

  //if the system if already small, dont create MG levels
  bool done = false;
  if(globalSize <= gCoarseSize){
    mg.AllocateLevelWorkSpace(mg.numLevels-1);
    coarse.setup(A, nullSpace, null, nullSpacePenalty);
    coarse.syncToDevice();
    mg.baseLevel = mg.numLevels-1;
    Lbase.syncToDevice();
    done = true;
  }

  //TODO: make the coarsen threasholds user-provided inputs
  // For now, let default to some sensible thresholds
  dfloat theta=0.0;
  if (mg.strtype==RUGESTUBEN) {
    theta=0.5; //default for 3D problems
    //See: A GPU accelerated aggregation algebraic multigrid method, R. Gandham, K. Esler, Y. Zhang.
  } else { // (type==SYMMETRIC)
    theta=0.08;
    //See: Algebraic Multigrid On Unstructured Meshes, P Vanek, J. Mandel, M. Brezina.
  }

  while(!done){
    /*Get current coarsest level*/
    amgLevel& L = mg.GetLevel<amgLevel>(mg.numLevels-1);

    /*Build smoother*/
    L.setupSmoother();

    /*Create new level*/
    amgLevel& Lcoarse = mg.AddLevel<amgLevel>();

    /* Coarsen level via AMG. Coarsen null vector */
    Lcoarse = coarsenAmgLevel(L, null,
                              mg.strtype, theta,
                              mg.aggtype);

    mg.AllocateLevelWorkSpace(mg.numLevels-2);
    L.syncToDevice();

    parCSR& Acoarse = Lcoarse.A;

    // Increase coarsening rate as we add levels.
    //See: Algebraic Multigrid On Unstructured Meshes, P Vanek, J. Mandel, M. Brezina.
    if (mg.strtype==SYMMETRIC)
      theta=theta/2;

    hlong globalCoarseSize;
    if (mg.coarsetype==COARSEEXACT) {
      globalCoarseSize = Acoarse.globalRowStarts[size];;
    } else { //COARSEOAS
      //OAS cares about Ncols for size
      globalCoarseSize = Acoarse.Ncols;
      Acoarse.comm.Allreduce(globalCoarseSize);
    }

    if(globalCoarseSize <= gCoarseSize || globalSize < 2*globalCoarseSize){
      LIBP_WARNING("AMG coarsening stalling, attemping coarse solver setup with dimension N=" << globalCoarseSize,
                   globalSize < 2*globalCoarseSize && rank==0);

      mg.AllocateLevelWorkSpace(mg.numLevels-1);
      Lcoarse.syncToDevice();
      coarse.setup(Acoarse, nullSpace, null, nullSpacePenalty);
      coarse.syncToDevice();
      mg.baseLevel = mg.numLevels-1;
      break;
    }
    globalSize = globalCoarseSize;
  }

  if(Comm::World().rank()==0) printf("done.\n");
}

} //namespace parAlmond

} //namespace libp
