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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsMatrix.hpp"
#include <random>
#include <algorithm>

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

namespace paradogs {

std::mt19937 RNG;

//------------------------------------------------------------------------
//
//  parCSR matrix
//
//------------------------------------------------------------------------

void parCSR::SpMV(const dfloat alpha, memory<dfloat>& x,
                  const dfloat beta, memory<dfloat>& y) {

  halo.ExchangeStart(x, 1);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  #pragma omp parallel for
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

  #pragma omp parallel for
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
  #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=diag.rowStarts[i]; jj<diag.rowStarts[i+1]; jj++)
      result += diag.vals[jj]*x[diag.cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }

  halo.ExchangeFinish(x, 1);

  #pragma omp parallel for
  for(dlong i=0; i<offd.nzRows; i++){ //local
    const dlong row = offd.rows[i];
    dfloat result = 0.0;
    for(dlong jj=offd.mRowStarts[i]; jj<offd.mRowStarts[i+1]; jj++)
      result += offd.vals[jj]*x[offd.cols[jj]];

    z[row] += alpha*result;
  }
}

//------------------------------------------------------------------------
//
//  parCSR matrix setup
//
//------------------------------------------------------------------------

//build a parCSR matrix from a distributed COO matrix
parCSR::parCSR(dlong _Nrows, dlong _Ncols,
               const dlong NNZ,
               memory<nonZero_t>& entries,
               const platform_t &_platform,
               comm_t _comm):
  platform(_platform),
  comm(_comm) {

  Nrows = _Nrows;
  Ncols = _Ncols;

  /*Get global row/col offsets*/
  hlong localNrows = static_cast<hlong>(Nrows);
  hlong localNcols = static_cast<hlong>(Ncols);
  comm.Scan(localNrows, rowOffsetU);
  comm.Scan(localNcols, colOffsetU);
  rowOffsetL = rowOffsetU-Nrows;
  colOffsetL = colOffsetU-Ncols;

  diag.rowStarts.malloc(Nrows+1);
  offd.rowStarts.malloc(Nrows+1);

  #pragma omp parallel for
  for (dlong n=0;n<Nrows+1;n++) {
    diag.rowStarts[n]=0;
    offd.rowStarts[n]=0;
  }

  // //count the entries in each row
  for (dlong n=0;n<NNZ;n++) {
    const dlong row = static_cast<dlong>(entries[n].row-rowOffsetL);
    if (   (entries[n].col <  colOffsetL)
        || (entries[n].col >= colOffsetU)) {
      offd.rowStarts[row+1]++;
    } else {
      diag.rowStarts[row+1]++;
    }
  }

  // count how many rows are shared
  offd.nzRows=0;
  for(dlong i=0; i<Nrows; i++)
    if (offd.rowStarts[i+1]>0) offd.nzRows++;

  offd.rows.malloc(offd.nzRows);
  offd.mRowStarts.malloc(offd.nzRows+1);

  // cumulative sum
  dlong cnt=0;
  offd.mRowStarts[0]=0;
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
  for (dlong n=0;n<NNZ;n++) {
    if (   (entries[n].col <  colOffsetL)
        || (entries[n].col >= colOffsetU)) {
      colIds[cnt++] = entries[n].col;
    }
  }
  haloSetup(colIds); //setup halo, and transform colIds to a local indexing

  // //fill the CSR matrices
  diag.cols.malloc(diag.nnz);
  offd.cols.malloc(offd.nnz);
  diag.vals.malloc(diag.nnz);
  offd.vals.malloc(offd.nnz);
  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for (dlong n=0;n<NNZ;n++) {
    if (   (entries[n].col <  colOffsetL)
        || (entries[n].col >= colOffsetU)) {
      offd.cols[offdCnt] = colIds[offdCnt];
      offd.vals[offdCnt] = entries[n].val;
      offdCnt++;
    } else {
      diag.cols[diagCnt] = static_cast<dlong>(entries[n].col-colOffsetL);
      diag.vals[diagCnt] = entries[n].val;
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


void parCSR::haloSetup(memory<hlong>& colIds) {

  //collect the unique nonlocal column ids
  memory<parallelId_t> parIds(offd.nnz);

  for (dlong n=0;n<offd.nnz;n++) {
    parIds[n].localId  = n;
    parIds[n].globalId = colIds[n];
  }

  //sort by global index
  sort(parIds.ptr(), parIds.ptr()+offd.nnz,
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
  sort(parIds.ptr(), parIds.ptr()+offd.nnz,
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
  for (dlong n=0; n<NlocalCols; n++)      colMap[n] = n+colOffsetL+1; //local rows
  for (dlong n=NlocalCols; n<Ncols; n++)  colMap[n] = -(offdcols[n-NlocalCols]+1);    //nonlocal rows

  //make a halo exchange to share column entries and an ogs for gsops accross columns
  bool verbose = false;
  halo.Setup(Ncols, colMap, comm, ogs::Pairwise, verbose, platform);

  //shift back to 0-indexed
  for (dlong n=0; n<Ncols; n++) colMap[n]=abs(colMap[n])-1;

  //update column numbering
  for (dlong n=0;n<offd.nnz;n++)
    colIds[n] = NlocalCols + parIds[n].newId;
}

//------------------------------------------------------------------------
//
//  parCSR Estimate max Eigenvalue of diagA^{-1}*A
//
//------------------------------------------------------------------------

dfloat parCSR::rhoDinvA(memory<dfloat>& null){

  int k = 10;

  hlong Ntotal = static_cast<hlong>(Nrows);
  comm.Allreduce(Ntotal);
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  memory<double> H(k*k, 0.0);

  // allocate memory for basis
  memory<dfloat> V((k+1)*Nrows);
  memory<dfloat> Vx(Ncols);

  /*Create rng*/
  std::uniform_real_distribution<dfloat> distrib(-0.5, 0.5);

  // generate a random vector for initial basis vector
  for(dlong n=0; n<Nrows; n++) Vx[n] = distrib(RNG);

  /*Project out null vector*/
  dfloat nulldot =0.0;
  for(dlong n=0; n<Nrows; n++) nulldot += null[n]*Vx[n];
  comm.Allreduce(nulldot);

  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++) Vx[n] -= nulldot*null[n];

  // dfloat norm_vo = vectorNorm(Nrows,Vx, comm);
  dfloat norm_vo=0.0;
  for(dlong n=0; n<Nrows; n++) norm_vo += Vx[n]*Vx[n];
  comm.Allreduce(norm_vo);
  norm_vo = sqrt(norm_vo);

  // vectorScale(Nrows, 1.0/norm_vo, Vx);
  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++) Vx[n] *= (1.0/norm_vo);

  //V[0] = Vx
  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++) V[n] = Vx[n];

  for(int j=0; j<k; j++){
    memory<dfloat> Vj   = V+j*Nrows;
    memory<dfloat> Vjp1 = V+(j+1)*Nrows;

    //Vx = V[j]
    #pragma omp parallel for
    for(dlong n=0; n<Nrows; n++) Vx[n] = Vj[n];

    // v[j+1] = invD*(A*v[j])
    SpMV(1.0, Vx, 0., Vjp1);

    // vectorDotStar(Nrows, diagInv, V[j+1]);
    #pragma omp parallel for
    for(dlong n=0; n<Nrows; n++) Vjp1[n] *= diagInv[n];

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      memory<dfloat> Vi = V+i*Nrows;
      // H(i,j) = v[i]'*A*v[j]
      // dfloat hij = vectorInnerProd(Nrows, V[i], V[j+1],comm);
      dfloat hij=0.0;
      for(dlong n=0; n<Nrows; n++) hij += Vi[n]*Vjp1[n];
      comm.Allreduce(hij);

      // v[j+1] = v[j+1] - hij*v[i]
      // vectorAdd(Nrows,-hij, V[i], 1.0, V[j+1]);
      #pragma omp parallel for
      for(dlong n=0; n<Nrows; n++) Vjp1[n] += -hij*Vi[n];

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){

      // dfloat norm_vj = vectorNorm(Nrows,V[j+1],comm);
      dfloat norm_vj=0.0;
      for(dlong n=0; n<Nrows; n++) norm_vj += Vjp1[n]*Vjp1[n];
      comm.Allreduce(norm_vj);
      norm_vj = sqrt(norm_vj);

      H[j+1+ j*k] = (double) norm_vj;

      // vectorScale(Nrows, 1./H[j+1 + j*k], V[j+1]);
      #pragma omp parallel for
      for(dlong n=0; n<Nrows; n++) Vjp1[n] *= (1./H[j+1 + j*k]);
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

} //namespace paradogs

} //namespace libp
