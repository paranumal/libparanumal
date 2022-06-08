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
#include "parAlmond/parAlmondparCSR.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace libp {

namespace parAlmond {

//------------------------------------------------------------------------
//
//  parCSR matrix
//
//------------------------------------------------------------------------

void parCSR::SpMV(const dfloat alpha, memory<dfloat>& x,
                  const dfloat beta, memory<dfloat>& y) {

  halo.ExchangeStart(x, 1);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=diag.rowStarts[i]; jj<diag.rowStarts[i+1]; jj++)
      result += diag.vals[jj]*x[diag.cols[jj]];

    if (beta!=0.0)
      y[i] = alpha*result + beta*y[i];
    else
      y[i] = alpha*result;
  }

  halo.ExchangeFinish(x, 1);

  // #pragma omp parallel for
  for(dlong i=0; i<offd.nzRows; i++){ //local
    const dlong row = offd.rows[i];
    dfloat result = 0.0;
    for(dlong jj=offd.mRowStarts[i]; jj<offd.mRowStarts[i+1]; jj++)
      result += offd.vals[jj]*x[offd.cols[jj]];

    y[row] += alpha*result;
  }
}

void parCSR::SpMV(const dfloat alpha, memory<dfloat>& x,
                  const dfloat beta, const memory<dfloat>& y, memory<dfloat>& z) {

  halo.ExchangeStart(x, 1);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=diag.rowStarts[i]; jj<diag.rowStarts[i+1]; jj++)
      result += diag.vals[jj]*x[diag.cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }

  halo.ExchangeFinish(x, 1);

  for(dlong i=0; i<offd.nzRows; i++){ //local
    const dlong row = offd.rows[i];
    dfloat result = 0.0;
    for(dlong jj=offd.mRowStarts[i]; jj<offd.mRowStarts[i+1]; jj++)
      result += offd.vals[jj]*x[offd.cols[jj]];

    z[row] += alpha*result;
  }
}

void parCSR::SpMV(const dfloat alpha, deviceMemory<dfloat>& o_x, const dfloat beta,
                  deviceMemory<dfloat>& o_y) {

  halo.ExchangeStart(o_x, 1);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (diag.NrowBlocks)
    SpMVcsrKernel1(diag.NrowBlocks, alpha, beta,
                   diag.o_blockRowStarts, diag.o_rowStarts,
                   diag.o_cols, diag.o_vals,
                   o_x, o_y);

  halo.ExchangeFinish(o_x, 1);

  const dfloat one = 1.0;
  if (offd.NrowBlocks)
    SpMVmcsrKernel(offd.NrowBlocks, alpha, one,
                   offd.o_blockRowStarts, offd.o_mRowStarts,
                   offd.o_rows, offd.o_cols, offd.o_vals,
                   o_x, o_y);
}

void parCSR::SpMV(const dfloat alpha, deviceMemory<dfloat>& o_x, const dfloat beta,
                  deviceMemory<dfloat>& o_y, deviceMemory<dfloat>& o_z) {

  halo.ExchangeStart(o_x, 1);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (diag.NrowBlocks)
    SpMVcsrKernel2(diag.NrowBlocks, alpha, beta,
                   diag.o_blockRowStarts, diag.o_rowStarts,
                   diag.o_cols, diag.o_vals,
                   o_x, o_y, o_z);

  halo.ExchangeFinish(o_x, 1);

  const dfloat one = 1.0;
  if (offd.NrowBlocks)
    SpMVmcsrKernel(offd.NrowBlocks, alpha, one,
                   offd.o_blockRowStarts, offd.o_mRowStarts,
                   offd.o_rows, offd.o_cols, offd.o_vals,
                   o_x, o_z);
}


//------------------------------------------------------------------------
//
//  parCSR matrix setup
//
//------------------------------------------------------------------------

//build a parCSR matrix from a distributed COO matrix
parCSR::parCSR(parCOO& A):       // number of nonzeros on this rank
  platform(A.platform),
  comm(A.comm) {

  int rank = comm.rank();
  // int size = comm.size();

  //copy global partition
  globalRowStarts = A.globalRowStarts;
  globalColStarts = A.globalColStarts;

  const hlong globalRowOffset = globalRowStarts[rank];
  const hlong globalColOffset = globalColStarts[rank];

  Nrows = static_cast<dlong>(globalRowStarts[rank+1]-globalRowStarts[rank]);
  Ncols = static_cast<dlong>(globalColStarts[rank+1]-globalColStarts[rank]);

  diag.rowStarts.malloc(Nrows+1, 0);
  offd.rowStarts.malloc(Nrows+1, 0);

  //count the entries in each row
  for (dlong n=0;n<A.nnz;n++) {
    const dlong row = (dlong) (A.entries[n].row - globalRowOffset);
    if (   (A.entries[n].col < globalColOffset)
        || (A.entries[n].col > globalColOffset+Ncols-1))
      offd.rowStarts[row+1]++;
    else
      diag.rowStarts[row+1]++;
  }

  offd.nzRows=0;

  // count how many rows are shared
  for(dlong i=0; i<Nrows; i++)
    if (offd.rowStarts[i+1]>0) offd.nzRows++;

  offd.rows.malloc(offd.nzRows);
  offd.mRowStarts.malloc(offd.nzRows+1);

  // cumulative sum
  dlong cnt=0;
  offd.mRowStarts[0] = 0;
  for(dlong i=0; i<Nrows; i++) {
    if (offd.rowStarts[i+1]>0) {
      offd.rows[cnt] = i; //record row id
      offd.mRowStarts[cnt+1] = offd.mRowStarts[cnt] + offd.rowStarts[i+1];
      cnt++;
    }
    diag.rowStarts[i+1] += diag.rowStarts[i];
    offd.rowStarts[i+1] += offd.rowStarts[i];
  }
  diag.nnz = diag.rowStarts[Nrows];
  offd.nnz = offd.rowStarts[Nrows];

  // Halo setup
  cnt=0;
  memory<hlong> colIds(offd.nnz);
  for (dlong n=0;n<A.nnz;n++) {
    if ( (A.entries[n].col < globalColOffset)
      || (A.entries[n].col > globalColOffset+Ncols-1))
      colIds[cnt++] = A.entries[n].col;
  }
  haloSetup(colIds); //setup halo, and transform colIds to a local indexing

  //fill the CSR matrices
  diag.cols.malloc(diag.nnz);
  offd.cols.malloc(offd.nnz);
  diag.vals.malloc(diag.nnz);
  offd.vals.malloc(offd.nnz);
  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for (dlong n=0;n<A.nnz;n++) {
    if ( (A.entries[n].col < globalColOffset)
      || (A.entries[n].col > globalColOffset+NlocalCols-1)) {
      offd.cols[offdCnt] = colIds[offdCnt];
      offd.vals[offdCnt] = A.entries[n].val;
      offdCnt++;
    } else {
      diag.cols[diagCnt] = static_cast<dlong>(A.entries[n].col - globalColOffset);
      diag.vals[diagCnt] = A.entries[n].val;
      diagCnt++;
    }
  }
}

//------------------------------------------------------------------------
//
//  parCSR halo setup
//
//------------------------------------------------------------------------

typedef struct {

  dlong localId;
  hlong globalId;

  dlong newId;

} parallelId_t;


void parCSR::haloSetup(memory<hlong> colIds) {

  int rank = comm.rank();

  const hlong globalOffset = globalColStarts[rank];

  //collect the unique nonlocal column ids
  memory<parallelId_t> parIds(offd.nnz);

  for (dlong n=0;n<offd.nnz;n++) {
    parIds[n].localId  = n;
    parIds[n].globalId = colIds[n];
  }

  //sort by global index
  std::sort(parIds.ptr(), parIds.ptr()+offd.nnz,
            [](const parallelId_t& a, const parallelId_t& b) {
              if(a.globalId < b.globalId) return true;
              if(a.globalId > b.globalId) return false;

              return (a.localId < b.localId);
            });

  //count unique nonlocal column ids
  dlong Noffdcols = 0; //number of unique columns
  if(offd.nnz) parIds[0].newId = Noffdcols;
  for (dlong n=1;n<offd.nnz;n++) {
    if (parIds[n].globalId != parIds[n-1].globalId)
      Noffdcols++;

    parIds[n].newId = Noffdcols;
  }
  if(offd.nnz) Noffdcols++;

  //record the global ids of the unique columns
  memory<hlong> offdcols(Noffdcols);
  Noffdcols = 0;
  if(offd.nnz) offdcols[Noffdcols++] = parIds[0].globalId;
  for (dlong n=1;n<offd.nnz;n++)
    if (parIds[n].globalId != parIds[n-1].globalId)
      offdcols[Noffdcols++] = parIds[n].globalId;

  //sort back to local order
  std::sort(parIds.ptr(), parIds.ptr()+offd.nnz,
            [](const parallelId_t& a, const parallelId_t& b) {
              if(a.localId < b.localId) return true;
              if(a.localId > b.localId) return false;

              return (a.globalId < b.globalId);
            });

  // be careful to make sure Ncols is set at this point
  NlocalCols = Ncols;
  Ncols += Noffdcols;

  //make an array of all the column ids required on this rank (local first)
  colMap.malloc(Ncols);
  for (dlong n=0; n<NlocalCols; n++)      colMap[n] = n+globalOffset+1; //local rows
  for (dlong n=NlocalCols; n<Ncols; n++)  colMap[n] = -(offdcols[n-NlocalCols]+1);    //nonlocal rows

  //make a halo exchange to share column entries and an ogs for gsops accross columns
  int verbose = 0;
  halo.Setup(Ncols, colMap, comm, ogs::Auto, verbose, platform);

  //shift back to 0-indexed
  for (dlong n=0; n<Ncols; n++) colMap[n]=std::abs(colMap[n])-1;

  //update column numbering
  for (dlong n=0;n<offd.nnz;n++)
    colIds[n] = NlocalCols + parIds[n].newId;
}

//------------------------------------------------------------------------
//
//  parCSR diagonal setup
//
//------------------------------------------------------------------------

void parCSR::diagSetup() {
  //fill the CSR matrices
  diagA.malloc(Ncols);
  diagInv.malloc(Ncols);

  for (dlong i=0;i<Nrows;i++) {
    const dlong start = diag.rowStarts[i];
    const dlong end   = diag.rowStarts[i+1];

    for (dlong j=start;j<end;j++) {
      //record the diagonal
      if (diag.cols[j]==i)
        diagA[i] = diag.vals[j];
    }
  }

  //fill the halo region
  halo.Exchange(diagA, 1);

  //compute the inverse diagonal
  for (dlong n=0;n<Nrows;n++)
    diagInv[n] = (diagA[n] != 0.0) ? 1.0/diagA[n] : 0.0;

  // estimate rho(invD * A)
  rho = rhoDinvA();
}


//------------------------------------------------------------------------
//
//  parCSR Estimate max Eigenvalue of diagA^{-1}*A
//
//------------------------------------------------------------------------

dfloat parCSR::rhoDinvA(){

  int size = comm.size();

  int k = 10;

  hlong Ntotal = globalRowStarts[size];
  if(k > Ntotal) k = static_cast<int>(Ntotal);

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  memory<double> H(k*k, 0.0);

  // allocate memory for basis
  memory<memory<dfloat>> V(k+1);
  memory<dfloat> Vx(Ncols);

  for(int i=0; i<=k; i++) {
    V[i].malloc(Nrows);
  }

  // generate a random vector for initial basis vector
  for(dlong n=0; n<Nrows; n++) Vx[n] = (dfloat) drand48();

  // dfloat norm_vo = vectorNorm(Nrows,Vx, comm);
  dfloat norm_vo=0.0;
  for(dlong n=0; n<Nrows; n++) norm_vo += Vx[n]*Vx[n];
  comm.Allreduce(norm_vo);
  norm_vo = sqrt(norm_vo);

  // vectorScale(Nrows, 1.0/norm_vo, Vx);
  for(dlong n=0; n<Nrows; n++) Vx[n] *= (1.0/norm_vo);

  //V[0] = Vx
  V[0].copyFrom(Vx, Nrows);

  for(int j=0; j<k; j++){
    //Vx = V[j]
    Vx.copyFrom(V[j], Nrows);

    // v[j+1] = invD*(A*v[j])
    SpMV(1.0, Vx, 0., V[j+1]);
    // vectorDotStar(Nrows, diagInv, V[j+1]);
    for(dlong n=0; n<Nrows; n++) V[j+1][n] *= diagInv[n];

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      // dfloat hij = vectorInnerProd(Nrows, V[i], V[j+1],comm);
      dfloat hij=0.0;
      for(dlong n=0; n<Nrows; n++) hij += V[i][n]*V[j+1][n];
      comm.Allreduce(hij);

      // v[j+1] = v[j+1] - hij*v[i]
      // vectorAdd(Nrows,-hij, V[i], 1.0, V[j+1]);
      for(dlong n=0; n<Nrows; n++) V[j+1][n] += -hij*V[i][n];

      H[i + j*k] = static_cast<double>(hij);
    }

    if(j+1 < k){

      // dfloat norm_vj = vectorNorm(Nrows,V[j+1],comm);
      dfloat norm_vj=0.0;
      for(dlong n=0; n<Nrows; n++) norm_vj += V[j+1][n]*V[j+1][n];
      comm.Allreduce(norm_vj);
      norm_vj = sqrt(norm_vj);

      H[j+1+ j*k] = static_cast<double>(norm_vj);

      // vectorScale(Nrows, 1./H[j+1 + j*k], V[j+1]);
      for(dlong n=0; n<Nrows; n++) V[j+1][n] *= (1./H[j+1 + j*k]);
    }
  }

  memory<double> WR(k);
  memory<double> WI(k);

  linAlg_t::matrixEigenValues(k, H, WR, WI);

  double RHO = 0.;

  for(int i=0; i<k; i++){
    double RHO_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(RHO < RHO_i) {
      RHO = RHO_i;
    }
  }

  // printf("weight = %g \n", RHO);

  return RHO;
}

void parCSR::syncToDevice() {

  if (Nrows) {
    //transfer matrix data
    diag.o_rowStarts = platform.malloc<dlong>(diag.rowStarts);

    diag.NrowBlocks=0;
    if (diag.nnz) {
      //setup row blocking
      dlong blockSum=0;
      diag.NrowBlocks=1;
      for (dlong i=0;i<Nrows;i++) {
        dlong rowSize = diag.rowStarts[i+1]-diag.rowStarts[i];

        //this may be pathalogically big. We can't currently run this
        LIBP_ABORT("Multiplicity of row: " << i << " is " << rowSize << " in parAlmond::parCSR setup and is too large.",
                   rowSize > parAlmond::NonzerosPerBlock);

        if (blockSum+rowSize > parAlmond::NonzerosPerBlock) { //adding this row will exceed the nnz per block
          diag.NrowBlocks++; //count the previous block
          blockSum=rowSize; //start a new row block
        } else {
          blockSum+=rowSize; //add this row to the block
        }
      }

      diag.blockRowStarts.malloc(diag.NrowBlocks+1, 0);

      blockSum=0;
      diag.NrowBlocks=1;
      for (dlong i=0;i<Nrows;i++) {
        dlong rowSize = diag.rowStarts[i+1]-diag.rowStarts[i];

        if (blockSum+rowSize > parAlmond::NonzerosPerBlock) { //adding this row will exceed the nnz per block
          diag.blockRowStarts[diag.NrowBlocks++] = i; //mark the previous block
          blockSum=rowSize; //start a new row block
        } else {
          blockSum+=rowSize; //add this row to the block
        }
      }
      diag.blockRowStarts[diag.NrowBlocks] = Nrows;
      diag.o_blockRowStarts = platform.malloc<dlong>(diag.blockRowStarts);

      //transfer matrix data
      diag.o_cols = platform.malloc<dlong>(diag.cols);
      diag.o_vals = platform.malloc<pfloat>(diag.vals);
    }

    if (offd.nnz) {
      //setup row blocking
      dlong blockSum=0;
      offd.NrowBlocks=1;
      for (dlong i=0;i<offd.nzRows;i++) {
        dlong rowSize = offd.mRowStarts[i+1]-offd.mRowStarts[i];

        //this row may be pathalogically big. We can't currently run this
        LIBP_ABORT("Multiplicity of row: " << i << " is " << rowSize << " in parAlmond::parCSR setup and is too large.",
                   rowSize > parAlmond::NonzerosPerBlock);

        if (blockSum+rowSize > parAlmond::NonzerosPerBlock) { //adding this row will exceed the nnz per block
          offd.NrowBlocks++; //count the previous block
          blockSum=rowSize; //start a new row block
        } else {
          blockSum+=rowSize; //add this row to the block
        }
      }

      offd.blockRowStarts.malloc(offd.NrowBlocks+1, 0);

      blockSum=0;
      offd.NrowBlocks=1;
      for (dlong i=0;i<offd.nzRows;i++) {
        dlong rowSize = offd.mRowStarts[i+1]-offd.mRowStarts[i];

        if (blockSum+rowSize > parAlmond::NonzerosPerBlock) { //adding this row will exceed the nnz per block
          offd.blockRowStarts[offd.NrowBlocks++] = i; //mark the previous block
          blockSum=rowSize; //start a new row block
        } else {
          blockSum+=rowSize; //add this row to the block
        }
      }
      offd.blockRowStarts[offd.NrowBlocks] = offd.nzRows;
      offd.o_blockRowStarts = platform.malloc<dlong>(offd.blockRowStarts);

      //transfer matrix data
      offd.o_rows       = platform.malloc<dlong>(offd.rows);
      offd.o_mRowStarts = platform.malloc<dlong>(offd.mRowStarts);

      offd.o_cols = platform.malloc<dlong>(offd.cols);
      offd.o_vals = platform.malloc<pfloat>(offd.vals);
    }

    if (diagA.size()) {
      o_diagA = platform.malloc<dfloat>(diagA);
      o_diagInv = platform.malloc<dfloat>(diagInv);
    }
  }
}

} //namespace parAlmond

} //namespace libp
