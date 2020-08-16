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

#ifndef PARALMOND_MATRIX_HPP
#define PARALMOND_MATRIX_HPP


namespace parAlmond {

struct nonZero_t {
  hlong row;
  hlong col;
  dfloat val;

  int ownerRank;
};

class matrix_t {

public:
  dlong Nrows;
  dlong Ncols;

  matrix_t(dlong N=0, dlong M=0): Nrows(N), Ncols(M) {};
  virtual ~matrix_t() {};

  virtual void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, dfloat *y)=0;
  virtual void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, const dfloat *y, dfloat *z)=0;
  virtual void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, const occa::memory o_y)=0;
  virtual void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, occa::memory o_y, occa::memory o_z)=0;
};

class CSR: public matrix_t {

public:
  dlong nnz;
  dlong  *rowStarts=NULL;
  dlong  *cols=NULL;
  dfloat *vals=NULL;

  occa::memory o_rowStarts;
  occa::memory o_cols;
  occa::memory o_vals;

  CSR(dlong N=0, dlong M=0);
  ~CSR();

  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, dfloat *y);
  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, const dfloat *y, dfloat *z);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, const occa::memory o_y);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, occa::memory o_y, occa::memory o_z);
};

class ELL: public matrix_t {

public:
  int   nnzPerRow;
  dlong actualNNZ;

  dlong  *cols=NULL;
  dfloat *vals=NULL;

  occa::memory o_cols;
  occa::memory o_vals;

  ELL(dlong N=0, dlong M=0);
  ~ELL();

  void syncToDevice(occa::device device);

  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, dfloat *y);
  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, const dfloat *y, dfloat *z);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, const occa::memory o_y);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, occa::memory o_y, occa::memory o_z);
};

class MCSR: public matrix_t {

public:
  dlong nnz;
  dlong actualRows;
  dlong  *rowStarts=NULL;
  dlong  *rows=NULL;
  dlong  *cols=NULL;
  dfloat *vals=NULL;

  occa::memory o_rowStarts;
  occa::memory o_rows;
  occa::memory o_cols;
  occa::memory o_vals;

  MCSR(dlong N=0, dlong M=0);
  ~MCSR();

  void syncToDevice(occa::device device);

  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, dfloat *y);
  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, const dfloat *y, dfloat *z);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, const occa::memory o_y);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, occa::memory o_y, occa::memory o_z);
};

class parCOO {

public:
  dlong nnz=0;
  nonZero_t *entries=NULL;
  hlong *globalStarts=NULL;

  parCOO(): nnz(0) {}
  ~parCOO() {
    if(entries) free(entries);
    if(globalStarts) free(globalStarts);
  }
};

class parCSR: public matrix_t {

public:
  platform_t& platform;

  //local
  CSR *diag;

  //non-local
  CSR *offd;

  dfloat *diagA=NULL;
  dfloat *diagInv=NULL;

  occa::memory o_diagA;
  occa::memory o_diagInv;

  bool nullSpace;
  dfloat nullSpacePenalty;
  dfloat *null=NULL;
  occa::memory o_null;

  //partition info
  MPI_Comm comm;
  hlong *globalRowStarts=NULL;
  hlong *globalColStarts=NULL;
  hlong *colMap=NULL;

  halo_t *halo = NULL;
  dlong NlocalCols = 0;

  parCSR(dlong N, dlong M, platform_t& _platform);
  parCSR(dlong N, dlong M, MPI_Comm Comm, platform_t& _platform);

  //build a parCSR matrix from a distributed COO matrix
  parCSR(dlong N,         // number of rows on this rank
         parCOO& A,
         bool NullSpace,          //switch for nullspace
         dfloat *Null,            //null vector (or low energy mode)
         dfloat NullSpacePenalty, //penalty parameter for rank boost
         MPI_Comm Comm,
         platform_t &_platform);

  ~parCSR();

  void haloSetup(hlong *colIds);

  dfloat rhoDinvA();

  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, dfloat *y);
  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, const dfloat *y, dfloat *z);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, const occa::memory o_y);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, occa::memory o_y, occa::memory o_z);
};


class parHYB: public matrix_t {

public:

  ELL  *E;
  MCSR *C;

  dfloat *diagA=NULL;
  dfloat *diagInv=NULL;

  occa::memory o_diagA;
  occa::memory o_diagInv;

  bool nullSpace;
  dfloat nullSpacePenalty;
  dfloat *null=NULL;
  occa::memory o_null;

  //partition info
  MPI_Comm comm;
  hlong *globalRowStarts=NULL;
  hlong *globalColStarts=NULL;
  hlong *colMap=NULL;

  halo_t *halo = NULL;
  dlong NlocalCols = 0;

  parHYB(dlong N=0, dlong M=0);
  parHYB(parCSR *A); //build from parCSR

  ~parHYB();

  void syncToDevice(platform_t& platform);

  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, dfloat *y);
  void SpMV(const dfloat alpha,        dfloat *x, const dfloat beta, const dfloat *y, dfloat *z);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, const occa::memory o_y);
  void SpMV(const dfloat alpha, occa::memory o_x, const dfloat beta, occa::memory o_y, occa::memory o_z);
};


} //namespace parAlmond

#endif