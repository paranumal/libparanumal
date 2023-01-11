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

#ifndef PARALMOND_PARCSR_HPP
#define PARALMOND_PARCSR_HPP

namespace libp {

namespace parAlmond {

class parCSR {
public:
  platform_t platform;
  comm_t comm;

  dlong Nrows=0;
  dlong Ncols=0;

  //local sparse matrix
  struct CSR {
    dlong nnz=0;
    dlong NrowBlocks=0;

    memory<dlong>  blockRowStarts;
    memory<dlong>  rowStarts;
    memory<dlong>  cols;
    memory<pfloat> vals;

    deviceMemory<dlong>  o_blockRowStarts;
    deviceMemory<dlong>  o_rowStarts;
    deviceMemory<dlong>  o_cols;
    deviceMemory<pfloat> o_vals;
  };
  CSR diag;

  //non-local sparse matrix
  struct MCSR {
    dlong nnz=0;
    dlong nzRows=0;
    dlong NrowBlocks=0;

    memory<dlong>  blockRowStarts;
    memory<dlong>  rowStarts;
    memory<dlong>  mRowStarts; //compressed version of rowStarts
    memory<dlong>  rows;
    memory<dlong>  cols;
    memory<pfloat> vals;

    deviceMemory<dlong>  o_blockRowStarts;
    deviceMemory<dlong>  o_mRowStarts;
    deviceMemory<dlong>  o_rows;
    deviceMemory<dlong>  o_cols;
    deviceMemory<pfloat> o_vals;
  };
  MCSR offd;

  memory<pfloat> diagA;
  memory<pfloat> diagInv;

  deviceMemory<pfloat> o_diagA;
  deviceMemory<pfloat> o_diagInv;

  //partition info
  memory<hlong> globalRowStarts;
  memory<hlong> globalColStarts;
  memory<hlong> colMap;

  ogs::halo_t halo;
  dlong NlocalCols = 0;

  //rho ~= cond(invD * A)
  pfloat rho=0.0;

  parCSR() = default;
  parCSR(dlong N, dlong M, platform_t& _platform, comm_t _comm):
    platform(_platform), comm(_comm), Nrows(N), Ncols(M) {}

  //build a parCSR matrix from a distributed COO matrix
  parCSR(parCOO& A);

  void haloSetup(memory<hlong> colIds);

  void diagSetup();

  pfloat rhoDinvA();

  void syncToDevice();

  void SpMV(const pfloat alpha, memory<pfloat>& x,
            const pfloat beta, memory<pfloat>& y);
  void SpMV(const pfloat alpha, memory<pfloat>& x,
            const pfloat beta, const memory<pfloat>& y, memory<pfloat>& z);

  void SpMV(const pfloat alpha, deviceMemory<pfloat>& o_x, const pfloat beta,
            deviceMemory<pfloat>& o_y);
  void SpMV(const pfloat alpha, deviceMemory<pfloat>& o_x, const pfloat beta,
            deviceMemory<pfloat>& o_y, deviceMemory<pfloat>& o_z);

  void smoothDampedJacobi(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_x,
                          const pfloat lambda, bool x_is_zero);

  void smoothChebyshev(deviceMemory<pfloat>& o_b, deviceMemory<pfloat>& o_x,
                       const pfloat lambda0, const pfloat lambda1,
                       bool x_is_zero, const int ChebyshevIterations);
};

} //namespace parAlmond

} //namespace libp

#endif
