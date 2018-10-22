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

#ifndef PARALMOND_EXACTSOLVE_HPP
#define PARALMOND_EXACTSOLVE_HPP

namespace parAlmond {

class exactSolver: public coarseSolver {

public:
  int coarseTotal;
  int coarseOffset;
  int *coarseOffsets=NULL;
  int *coarseCounts=NULL;

  int N;
  dfloat *invCoarseA=NULL;
  occa::memory o_invCoarseAT;

  dfloat *xLocal=NULL;
  dfloat *rhsLocal=NULL;

  dfloat *xCoarse=NULL;
  dfloat *rhsCoarse=NULL;
  occa::memory o_rhsCoarse;

  exactSolver(setupAide options, MPI_Comm comm_);
  ~exactSolver();

  int getTargetSize();

  void setup(parCSR *A);

  void Report(int lev);

  void solve(dfloat *rhs, dfloat *x);
  void solve(occa::memory o_rhs, occa::memory o_x);
};

}

#endif