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

#ifndef PARALMOND_PARCSR_HPP
#define PARALMOND_PARCSR_HPP

namespace parAlmond {

class parCSR {
public:
  platform_t& platform;
  MPI_Comm comm;

  dlong Nrows;
  dlong Ncols;

  //local sparse matrix
  struct CSR {
    dlong nnz=0;
    dlong NrowBlocks=0;

    dlong  *blockRowStarts=nullptr;
    dlong  *rowStarts=nullptr;
    dlong  *cols=nullptr;
    pfloat *vals=nullptr;

    occa::memory o_blockRowStarts;
    occa::memory o_rowStarts;
    occa::memory o_cols;
    occa::memory o_vals;
  };
  CSR diag;

  //non-local sparse matrix
  struct MCSR {
    dlong nnz=0;
    dlong nzRows=0;
    dlong NrowBlocks=0;

    dlong  *blockRowStarts=nullptr;
    dlong  *rowStarts=nullptr;
    dlong  *mRowStarts=nullptr; //compressed version of rowStarts
    dlong  *rows=nullptr;
    dlong  *cols=nullptr;
    pfloat *vals=nullptr;

    occa::memory o_blockRowStarts;
    occa::memory o_mRowStarts;
    occa::memory o_rows;
    occa::memory o_cols;
    occa::memory o_vals;
  };
  MCSR offd;

  pfloat *diagA=nullptr;
  pfloat *diagInv=nullptr;

  occa::memory o_diagA;
  occa::memory o_diagInv;

  //partition info
  hlong *globalRowStarts=nullptr;
  hlong *globalColStarts=nullptr;
  hlong *colMap=nullptr;

  halo_t *halo = nullptr;
  dlong NlocalCols = 0;

  //rho ~= cond(invD * A)
  pfloat rho=0.0;

  parCSR(dlong N, dlong M, platform_t& _platform, MPI_Comm _comm):
    platform(_platform), comm(_comm), Nrows(N), Ncols(M) {}

  //build a parCSR matrix from a distributed COO matrix
  parCSR(parCOO& A);

  ~parCSR();

  void haloSetup(hlong *colIds);

  void diagSetup();

  pfloat rhoDinvA();

  void syncToDevice();

  void SpMV(const pfloat alpha, pfloat *x,
            const pfloat beta, pfloat *y);
  void SpMV(const pfloat alpha, pfloat *x,
            const pfloat beta, const pfloat *y, pfloat *z);

  void SpMV(const pfloat alpha, occa::memory& o_x, const pfloat beta,
            occa::memory& o_y);
  void SpMV(const pfloat alpha, occa::memory& o_x, const pfloat beta,
            occa::memory& o_y, occa::memory& o_z);

  void smoothDampedJacobi(occa::memory& o_r, occa::memory& o_x,
                          const pfloat lambda, bool x_is_zero,
                          occa::memory& o_scratch);

  void smoothChebyshev(occa::memory& o_b, occa::memory& o_x,
                       const pfloat lambda0, const pfloat lambda1,
                       bool x_is_zero, occa::memory& o_scratch,
                       const int ChebyshevIterations);
};



} //namespace parAlmond

#endif
