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
#include "parAlmond/parAlmondparCSR.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace parAlmond {

//------------------------------------------------------------------------
//
//  parCSR matrix
//
//------------------------------------------------------------------------

void parCSR::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, dfloat *y) {

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

  halo->Exchange(x, 1, ogs_dfloat);

  // #pragma omp parallel for
  for(dlong i=0; i<offd.nzRows; i++){ //local
    const dlong row = offd.rows[i];
    dfloat result = 0.0;
    for(dlong jj=offd.rowStarts[i]; jj<offd.rowStarts[i+1]; jj++)
      result += offd.vals[jj]*x[offd.cols[jj]];

    y[row] += alpha*result;
  }
}

void parCSR::SpMV(const dfloat alpha, dfloat *x,
                  const dfloat beta, const dfloat *y, dfloat *z) {

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){ //local
    dfloat result = 0.0;
    for(dlong jj=diag.rowStarts[i]; jj<diag.rowStarts[i+1]; jj++)
      result += diag.vals[jj]*x[diag.cols[jj]];

    z[i] = alpha*result + beta*y[i];
  }

  halo->Exchange(x, 1, ogs_dfloat);

  for(dlong i=0; i<offd.nzRows; i++){ //local
    const dlong row = offd.rows[i];
    dfloat result = 0.0;
    for(dlong jj=offd.rowStarts[i]; jj<offd.rowStarts[i+1]; jj++)
      result += offd.vals[jj]*x[offd.cols[jj]];

    z[row] += alpha*result;
  }
}

void parCSR::SpMV(const dfloat alpha, occa::memory& o_x, const dfloat beta,
                  occa::memory& o_y) {

  halo->ExchangeStart(o_x, 1, ogs_dfloat);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (Nrows)
    SpMVcsrKernel1(Nrows, alpha, beta,
                   diag.o_rowStarts, diag.o_cols, diag.o_vals,
                   o_x, o_y);

  halo->ExchangeFinish(o_x, 1, ogs_dfloat);

  const dfloat one = 1.0;
  if (offd.nzRows)
    SpMVmcsrKernel(offd.nzRows, alpha, one,
                    offd.o_rowStarts, offd.o_rows, offd.o_cols, offd.o_vals,
                    o_x, o_y);
}

void parCSR::SpMV(const dfloat alpha, occa::memory& o_x, const dfloat beta,
                  occa::memory& o_y, occa::memory& o_z) {

  halo->ExchangeStart(o_x, 1, ogs_dfloat);

  // z[i] = beta*y[i] + alpha* (sum_{ij} Aij*x[j])
  if (Nrows)
    SpMVcsrKernel2(Nrows, alpha, beta,
                   diag.o_rowStarts, diag.o_cols, diag.o_vals,
                   o_x, o_y, o_z);

  halo->ExchangeFinish(o_x, 1, ogs_dfloat);

  const dfloat one = 1.0;
  if (offd.nzRows)
    SpMVmcsrKernel(offd.nzRows, alpha, one,
                    offd.o_rowStarts, offd.o_rows, offd.o_cols, offd.o_vals,
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

  int rank;
  int size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  Nrows = (dlong)(A.globalStarts[rank+1]-A.globalStarts[rank]);
  Ncols = Nrows;

  //copy global partition
  globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  memcpy(globalRowStarts, A.globalStarts, (size+1)*sizeof(hlong));
  memcpy(globalColStarts, A.globalStarts, (size+1)*sizeof(hlong));

  const hlong globalOffset = globalRowStarts[rank];

  diag.rowStarts = (dlong *) calloc(Nrows+1, sizeof(dlong));

  int* offdRowCounts = (dlong *) calloc(Nrows+1, sizeof(dlong));

  //count the entries in each row
  for (dlong n=0;n<A.nnz;n++) {
    const dlong row = (dlong) (A.entries[n].row - globalOffset);
    if (   (A.entries[n].col < globalOffset)
        || (A.entries[n].col > globalOffset+Nrows-1))
      offdRowCounts[row+1]++;
    else
      diag.rowStarts[row+1]++;
  }

  offd.nzRows=0;

  // count how many rows are shared
  for(dlong i=0; i<Nrows; i++)
    if (offdRowCounts[i+1]>0) offd.nzRows++;

  offd.rows      = (dlong *) calloc(offd.nzRows, sizeof(dlong));
  offd.rowStarts = (dlong *) calloc(offd.nzRows+1, sizeof(dlong));

  // cumulative sum
  dlong cnt=0;
  for(dlong i=0; i<Nrows; i++) {

    diag.rowStarts[i+1] += diag.rowStarts[i];

    if (offdRowCounts[i+1]>0) {
      offd.rows[cnt] = i; //record row id
      offd.rowStarts[cnt+1] = offd.rowStarts[cnt] + offdRowCounts[i+1];
      cnt++;
    }
  }
  diag.nnz = diag.rowStarts[Nrows];
  offd.nnz = offd.rowStarts[offd.nzRows];

  free(offdRowCounts);

  // Halo setup
  cnt=0;
  hlong *colIds = (hlong *) malloc(offd.nnz*sizeof(hlong));
  for (dlong n=0;n<A.nnz;n++) {
    if ( (A.entries[n].col < globalOffset)
      || (A.entries[n].col > globalOffset+Nrows-1))
      colIds[cnt++] = A.entries[n].col;
  }
  haloSetup(colIds); //setup halo, and transform colIds to a local indexing

  //fill the CSR matrices
  diagA   = (dfloat *) calloc(Ncols, sizeof(dfloat));
  diagInv = (dfloat *) calloc(Ncols, sizeof(dfloat));
  diag.cols = (dlong *)  calloc(diag.nnz, sizeof(dlong));
  offd.cols = (dlong *)  calloc(offd.nnz, sizeof(dlong));
  diag.vals = (dfloat *) calloc(diag.nnz, sizeof(dfloat));
  offd.vals = (dfloat *) calloc(offd.nnz, sizeof(dfloat));
  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for (dlong n=0;n<A.nnz;n++) {
    if ( (A.entries[n].col < globalOffset)
      || (A.entries[n].col > globalOffset+Nrows-1)) {
      offd.cols[offdCnt] = colIds[offdCnt];
      offd.vals[offdCnt] = A.entries[n].val;
      offdCnt++;
    } else {
      diag.cols[diagCnt] = (dlong) (A.entries[n].col - globalOffset);
      diag.vals[diagCnt] = A.entries[n].val;

      //record the diagonal
      dlong row = (dlong) (A.entries[n].row - globalOffset);
      if (row==diag.cols[diagCnt])
        diagA[row] = diag.vals[diagCnt];

      diagCnt++;
    }
  }
  free(colIds);

  //fill the halo region
  halo->Exchange(diagA, 1, ogs_dfloat);

  //compute the inverse diagonal
  for (dlong n=0;n<Nrows;n++) diagInv[n] = 1.0/diagA[n];
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

// compare on global indices
static int CompareGlobalId(const void *a, const void *b){

  parallelId_t *fa = (parallelId_t*) a;
  parallelId_t *fb = (parallelId_t*) b;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  return 0;
}

// compare on local indices
static int CompareLocalId(const void *a, const void *b){

  parallelId_t *fa = (parallelId_t*) a;
  parallelId_t *fb = (parallelId_t*) b;

  if(fa->localId < fb->localId) return -1;
  if(fa->localId > fb->localId) return +1;

  if(fa->globalId < fb->globalId) return -1;
  if(fa->globalId > fb->globalId) return +1;

  return 0;
}

void parCSR::haloSetup(hlong *colIds) {

  int rank;
  MPI_Comm_rank(comm, &rank);

  const hlong globalOffset = globalColStarts[rank];

  //collect the unique nonlocal column ids
  parallelId_t*  parIds = (parallelId_t*) malloc(offd.nnz*sizeof(parallelId_t));

  for (dlong n=0;n<offd.nnz;n++) {
    parIds[n].localId  = n;
    parIds[n].globalId = colIds[n];
  }

  //sort by global index
  qsort(parIds, offd.nnz, sizeof(parallelId_t), CompareGlobalId);

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
  hlong *offdcols = (hlong *) malloc(Noffdcols*sizeof(hlong));
  Noffdcols = 0;
  if(offd.nnz) offdcols[Noffdcols++] = parIds[0].globalId;
  for (dlong n=1;n<offd.nnz;n++)
    if (parIds[n].globalId != parIds[n-1].globalId)
      offdcols[Noffdcols++] = parIds[n].globalId;

  //sort back to local order
  qsort(parIds, offd.nnz, sizeof(parallelId_t), CompareLocalId);

  // be careful to make sure Ncols is set at this point
  NlocalCols = Ncols;
  Ncols += Noffdcols;

  //make an array of all the column ids required on this rank (local first)
  colMap = (hlong*) malloc(Ncols*sizeof(hlong));
  for (dlong n=0; n<NlocalCols; n++)      colMap[n] = n+globalOffset+1; //local rows
  for (dlong n=NlocalCols; n<Ncols; n++)  colMap[n] = -(offdcols[n-NlocalCols]+1);    //nonlocal rows

  //make a halo exchange to share column entries and an ogs for gsops accross columns
  int verbose = 0;
  halo = halo_t::Setup(Ncols, colMap, comm, verbose, platform);

  //shift back to 0-indexed
  for (dlong n=0; n<Ncols; n++) colMap[n]=abs(colMap[n])-1;

  //update column numbering
  for (dlong n=0;n<offd.nnz;n++)
    colIds[n] = NlocalCols + parIds[n].newId;

  free(parIds);
}

parCSR::~parCSR() {
  if (diag.blockRowStarts) free(diag.blockRowStarts);
  if (diag.rowStarts) free(diag.rowStarts);
  if (diag.cols) free(diag.cols);
  if (diag.vals) free(diag.vals);

  if (offd.blockRowStarts) free(offd.blockRowStarts);
  if (offd.rowStarts) free(offd.rowStarts);
  if (offd.rows) free(offd.rows);
  if (offd.cols) free(offd.cols);
  if (offd.vals) free(offd.vals);

  if (diagA) free(diagA);
  if (diagInv) free(diagInv);

  if (o_diagA.size()) o_diagA.free();
  if (o_diagInv.size()) o_diagInv.free();

  if (globalRowStarts) free(globalRowStarts);
  if (globalColStarts) free(globalColStarts);
  if (colMap) free(colMap);

  if (halo)   halo->Free();
}

//------------------------------------------------------------------------
//
//  parCSR Estimate max Eigenvalue of diagA^{-1}*A
//
//------------------------------------------------------------------------

dfloat parCSR::rhoDinvA(){

  int size;
  MPI_Comm_size(comm, &size);

  int k = 10;

  hlong Ntotal = globalRowStarts[size];
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double *H = (double *) calloc(k*k,sizeof(double));

  // allocate memory for basis
  dfloat **V = (dfloat **) calloc(k+1, sizeof(dfloat *));
  dfloat *Vx = (dfloat *) calloc(Ncols, sizeof(dfloat));

  for(int i=0; i<=k; i++)
    V[i] = (dfloat *) calloc(Nrows, sizeof(dfloat));

  // generate a random vector for initial basis vector
  for(dlong n=0; n<Nrows; n++) Vx[n] = (dfloat) drand48();

  // dfloat norm_vo = vectorNorm(Nrows,Vx, comm);
  dfloat norm_vo=0.0, gnorm_vo=0.0;
  for(dlong n=0; n<Nrows; n++) norm_vo += Vx[n]*Vx[n];
  MPI_Allreduce(&norm_vo, &gnorm_vo, 1, MPI_DFLOAT, MPI_SUM, comm);
  norm_vo = sqrt(gnorm_vo);

  // vectorScale(Nrows, 1.0/norm_vo, Vx);
  for(dlong n=0; n<Nrows; n++) Vx[n] *= (1.0/norm_vo);

  //V[0] = Vx
  memcpy(V[0], Vx, Nrows*sizeof(dfloat));

  for(int j=0; j<k; j++){
    //Vx = V[j]
    memcpy(Vx, V[j], Nrows*sizeof(dfloat));

    // v[j+1] = invD*(A*v[j])
    SpMV(1.0, Vx, 0., V[j+1]);
    // vectorDotStar(Nrows, diagInv, V[j+1]);
    for(dlong n=0; n<Nrows; n++) V[j+1][n] *= diagInv[n];

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      // dfloat hij = vectorInnerProd(Nrows, V[i], V[j+1],comm);
      dfloat local_hij=0.0, hij=0.0;
      for(dlong n=0; n<Nrows; n++) local_hij += V[i][n]*V[j+1][n];
      MPI_Allreduce(&local_hij, &hij, 1, MPI_DFLOAT, MPI_SUM, comm);

      // v[j+1] = v[j+1] - hij*v[i]
      // vectorAdd(Nrows,-hij, V[i], 1.0, V[j+1]);
      for(dlong n=0; n<Nrows; n++) V[j+1][n] += -hij*V[i][n];

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){

      // dfloat norm_vj = vectorNorm(Nrows,V[j+1],comm);
      dfloat norm_vj=0.0, gnorm_vj=0.0;
      for(dlong n=0; n<Nrows; n++) norm_vj += V[j+1][n]*V[j+1][n];
      MPI_Allreduce(&norm_vj, &gnorm_vj, 1, MPI_DFLOAT, MPI_SUM, comm);
      norm_vj = sqrt(gnorm_vj);

      H[j+1+ j*k] = (double) norm_vj;

      // vectorScale(Nrows, 1./H[j+1 + j*k], V[j+1]);
      for(dlong n=0; n<Nrows; n++) V[j+1][n] *= (1./H[j+1 + j*k]);
    }
  }

  double *WR = (double *) calloc(k,sizeof(double));
  double *WI = (double *) calloc(k,sizeof(double));

  matrixEigenValues(k, H, WR, WI);

  double rho = 0.;

  for(int i=0; i<k; i++){
    double rho_i  = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);

    if(rho < rho_i) {
      rho = rho_i;
    }
  }

  free(H);
  free(WR);
  free(WI);

  // free memory
  for(int i=0; i<=k; i++) free(V[i]);
  free(Vx);
  free(V);

  // printf("weight = %g \n", rho);

  return rho;
}

void parCSR::syncToDevice() {

  if (Nrows) {
    diag.o_rowStarts = platform.malloc((Nrows+1)*sizeof(dlong), diag.rowStarts);

    if (diag.nnz) {
      diag.o_cols = platform.malloc(diag.nnz*sizeof(dlong),   diag.cols);
      diag.o_vals = platform.malloc(diag.nnz*sizeof(dfloat),  diag.vals);
    }

    o_diagA   = platform.malloc(Nrows*sizeof(dfloat), diagA);
    o_diagInv = platform.malloc(Nrows*sizeof(dfloat), diagInv);

    if (offd.nzRows) {
      offd.o_rows      = platform.malloc(offd.nzRows*sizeof(dlong), offd.rows);
      offd.o_rowStarts = platform.malloc((offd.nzRows+1)*sizeof(dlong), offd.rowStarts);
    }
    if (offd.nnz) {
      offd.o_cols = platform.malloc(offd.nnz*sizeof(dlong),   offd.cols);
      offd.o_vals = platform.malloc(offd.nnz*sizeof(dfloat),  offd.vals);
    }
  }
}

} //namespace parAlmond
