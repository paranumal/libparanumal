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

#ifndef PARADOGS_MATRIX_HPP
#define PARADOGS_MATRIX_HPP 1

#include "parAdogs.hpp"
#include "ogs.hpp"

namespace libp {

namespace paradogs {

struct nonZero_t {
  hlong row;
  hlong col;
  dfloat val;
};


class parCSR {
public:
  platform_t platform;
  comm_t comm;

  dlong Nrows=0;
  dlong Ncols=0;

  //partition info
  hlong rowOffsetL=0, rowOffsetU=0;
  hlong colOffsetL=0, colOffsetU=0;

  //local sparse matrix
  struct CSR {
    dlong nnz=0;
    memory<dlong>  rowStarts;
    memory<dlong>  cols;
    memory<pfloat> vals;
  };
  CSR diag;

  //non-local sparse matrix
  struct MCSR {
    dlong nnz=0;
    dlong nzRows=0;

    memory<dlong>  rowStarts;
    memory<dlong>  mRowStarts;
    memory<dlong>  rows;
    memory<dlong>  cols;
    memory<pfloat> vals;
  };
  MCSR offd;

  memory<dfloat> diagA;
  memory<dfloat> diagInv;

  /*communcation info*/
  dlong NlocalCols = 0;
  ogs::halo_t halo;
  memory<hlong> colMap;

  //rho ~= cond(invD * A)
  dfloat rho=0.0;

  parCSR()=default;
  parCSR(dlong N, dlong M, platform_t& _platform, comm_t _comm):
    platform(_platform), comm(_comm), Nrows(N), Ncols(M) {}

  //build a parCSR matrix from a distributed COO matrix
  parCSR(dlong _Nrows, dlong _Ncols,
         const dlong NNZ,
         memory<nonZero_t>& entries,
         const platform_t &_platform,
         comm_t comm);

  void haloSetup(memory<hlong>& colIds);

  // estimate rho(invD * A)
  dfloat rhoDinvA(memory<dfloat>& null);

  /*Aggregate via distance-2 PMIS*/
  void Aggregate(dlong& cNverts,
                 const dfloat theta,
                 memory<hlong>& FineToCoarse);

  void GalerkinProduct(const parCSR &A, const parCSR &P);

  void SpMV(const dfloat alpha, memory<dfloat>& x,
            const dfloat beta, memory<dfloat>& y);
  void SpMV(const dfloat alpha, memory<dfloat>& x,
            const dfloat beta, const memory<dfloat>& y, memory<dfloat>& z);

  void SmoothChebyshev(memory<dfloat>& b, memory<dfloat>& x,
                       const dfloat lambda0, const dfloat lambda1,
                       const bool xIsZero, memory<dfloat>& scratch,
                       const int ChebyshevIterations);
};

} //namespace paradogs

} //namespace libp

#endif

