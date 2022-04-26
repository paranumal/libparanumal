/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifndef PARADOGS_MULTIGRID_HPP
#define PARADOGS_MULTIGRID_HPP 1

#include "parAdogs.hpp"
#include "parAdogs/parAdogsMatrix.hpp"

namespace libp {

namespace paradogs {

class mgLevel_t {
public:
  dlong Nrows=0, Ncols=0;
  hlong Nglobal=0;

  parCSR A, P, R;

  /*null vector*/
  memory<dfloat> null;

  /*Fiedler vector*/
  memory<dfloat> Fiedler;

  /*Vcycle storage*/
  memory<dfloat> RHS;
  memory<dfloat> X;
  memory<dfloat> RES;
  memory<dfloat> scratch;

  dfloat lambda1, lambda0; //smoothing params

  /*Create graph Laplacian*/
  void CreateLaplacian(const dlong Nelements,
                       const int Nfaces,
                       const memory<hlong>& EToE,
                       comm_t comm);

  /*Construct a coarse level*/
  void CoarsenLevel(mgLevel_t &Lf, const dfloat theta);

  void SetupSmoother();

  void AllocateScratch(const int l);

  /*Compute Fiedler vector directly*/
  void FiedlerVector();

  /*Multigrid functions*/
  void Smooth(memory<dfloat>& r, memory<dfloat>& x, const bool xIsZero);
  void Residual(memory<dfloat>& r, memory<dfloat>& x, memory<dfloat>& res);
  void Coarsen(memory<dfloat>& x, memory<dfloat>& xC);
  void Prolongate(memory<dfloat>& xC, memory<dfloat>& x);
};

parCSR TentativeProlongator(const dlong Nf,
                            const dlong Nc,
                            platform_t& platform,
                            comm_t comm,
                            memory<hlong>& FineToCoarse,
                            memory<dfloat>& FineNull,
                            memory<dfloat>& CoarseNull);

parCSR SmoothProlongator(const parCSR& A,
                         const parCSR& T);

parCSR Transpose(const parCSR& A);

parCSR SpMM(const parCSR& A, const parCSR& B);

class coarseSolver_t {

public:
  comm_t comm;

  int N=0;
  int Nrows=0;
  int Ncols=0;

  int coarseTotal=0;
  memory<int> coarseCounts;
  memory<int> coarseOffsets;

  memory<dfloat> invA;
  memory<dfloat> grhs;

  void Setup(parCSR& A, memory<dfloat>& null);
  void Solve(memory<dfloat>& r, memory<dfloat>& x);
};

} //namespace paradogs

} //namespace libp

#endif

