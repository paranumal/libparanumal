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

namespace parAlmond {

//------------------------------------------------------------------------
//
//  CSR matrix
//
//------------------------------------------------------------------------

CSR::CSR(dlong N, dlong M): matrix_t(N,M) {}

CSR::~CSR() {
  free(rowStarts);
  free(cols);
  free(vals);

  if (o_rowStarts.size()) o_rowStarts.free();
  if (o_cols.size()) o_cols.free();
  if (o_vals.size()) o_vals.free();
}

void CSR::syncToDevice(occa::device& device) {
  if (Nrows) {
    o_rowStarts = device.malloc((Nrows+1)*sizeof(dlong), rowStarts);
  }
  if (nnz) {
    o_cols = device.malloc(nnz*sizeof(dlong),   cols);
    o_vals = device.malloc(nnz*sizeof(dfloat),  vals);
  }
}

//------------------------------------------------------------------------
//
//  MCSR matrix
//
//------------------------------------------------------------------------
MCSR::MCSR(dlong N, dlong M): matrix_t(N,M) {}

MCSR::~MCSR() {
  free(rowStarts);
  free(rows);
  free(cols);
  free(vals);

  if (o_rowStarts.size()) o_rowStarts.free();
  if (o_rows.size()) o_rows.free();
  if (o_cols.size()) o_cols.free();
  if (o_vals.size()) o_vals.free();
}

void MCSR::syncToDevice(occa::device& device) {
  if (actualRows) {
    o_rowStarts = device.malloc((actualRows+1)*sizeof(dlong), rowStarts);
    o_rows      = device.malloc(actualRows*sizeof(dlong), rows);
  }
  if (nnz) {
    o_cols = device.malloc(nnz*sizeof(dlong),   cols);
    o_vals = device.malloc(nnz*sizeof(dfloat),  vals);
  }
}

//------------------------------------------------------------------------
//
//  parCSR matrix
//
//------------------------------------------------------------------------
parCSR::parCSR(dlong N, dlong M, platform_t& _platform):
  matrix_t(N,M), platform(_platform) {
  diag = new CSR(N,M);
  offd = new CSR(N,M);

  nullSpace=false;
}

parCSR::parCSR(dlong N, dlong M,
               MPI_Comm comm_,
               platform_t& _platform): matrix_t(N,M), platform(_platform) {
  MPI_Comm_dup(comm_, &comm);

  diag = new CSR(N,M);
  offd = new CSR(N,M);

  nullSpace=false;
}

//build a parCSR matrix from a distributed COO matrix (assumes square)
parCSR::parCSR(dlong N,         // number of rows on this rank
               parCOO& A,       // number of nonzeros on this rank
               bool NullSpace,          //switch for nullspace
               dfloat *Null,            //null vector (or low energy mode)
               dfloat NullSpacePenalty, //penalty parameter for rank boost
               MPI_Comm comm_,
               platform_t& _platform): platform(_platform) {

  MPI_Comm_dup(comm_, &comm);

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  Nrows = N;
  Ncols = N;

  //copy global partition
  globalRowStarts = (hlong *) calloc(size+1,sizeof(hlong));
  globalColStarts = (hlong *) calloc(size+1,sizeof(hlong));
  memcpy(globalRowStarts, A.globalStarts, (size+1)*sizeof(hlong));
  memcpy(globalColStarts, A.globalStarts, (size+1)*sizeof(hlong));

  hlong globalOffset = globalRowStarts[rank];

  null = (dfloat *) calloc(Nrows, sizeof(dfloat));
  memcpy(null, Null, Nrows*sizeof(dfloat));

  nullSpace = NullSpace;
  nullSpacePenalty = NullSpacePenalty;

  diag = new CSR(Nrows,Nrows);
  offd = new CSR(Nrows,Nrows);

  diag->rowStarts = (dlong *) calloc(Nrows+1, sizeof(dlong));
  offd->rowStarts = (dlong *) calloc(Nrows+1, sizeof(dlong));

  //count the entries in each row
  for (dlong n=0;n<A.nnz;n++) {
    dlong row = (dlong) (A.entries[n].row - globalOffset);
    if (   (A.entries[n].col < globalOffset)
        || (A.entries[n].col > globalOffset+Nrows-1))
      offd->rowStarts[row+1]++;
    else
      diag->rowStarts[row+1]++;
  }

  // cumulative sum
  for(dlong i=0; i<Nrows; i++) {
    diag->rowStarts[i+1] += diag->rowStarts[i];
    offd->rowStarts[i+1] += offd->rowStarts[i];
  }
  diag->nnz = diag->rowStarts[Nrows];
  offd->nnz = offd->rowStarts[Nrows];

  // Halo setup
  hlong *colIds = (hlong *) malloc(offd->nnz*sizeof(hlong));
  dlong cnt=0;
  for (dlong n=0;n<A.nnz;n++) {
    if ( (A.entries[n].col < globalOffset)
      || (A.entries[n].col > globalOffset+N-1))
      colIds[cnt++] = A.entries[n].col;
  }
  this->haloSetup(colIds);

  //fill the CSR matrices
  diagA   = (dfloat *) calloc(Ncols, sizeof(dfloat));
  diagInv = (dfloat *) calloc(Ncols, sizeof(dfloat));
  diag->cols = (dlong *)  calloc(diag->nnz, sizeof(dlong));
  offd->cols = (dlong *)  calloc(offd->nnz, sizeof(dlong));
  diag->vals = (dfloat *) calloc(diag->nnz, sizeof(dfloat));
  offd->vals = (dfloat *) calloc(offd->nnz, sizeof(dfloat));
  dlong diagCnt = 0;
  dlong offdCnt = 0;
  for (dlong n=0;n<A.nnz;n++) {
    if ( (A.entries[n].col < globalOffset)
      || (A.entries[n].col > globalOffset+Nrows-1)) {
      offd->cols[offdCnt] = colIds[offdCnt];
      offd->vals[offdCnt] = A.entries[n].val;
      offdCnt++;
    } else {
      diag->cols[diagCnt] = (dlong) (A.entries[n].col - globalOffset);
      diag->vals[diagCnt] = A.entries[n].val;

      //record the diagonal
      dlong row = (dlong) (A.entries[n].row - globalOffset);
      if (row==diag->cols[diagCnt])
        diagA[row] = diag->vals[diagCnt];

      diagCnt++;
    }
  }

  //fill the halo region
  halo->Exchange(diagA, 1, ogs_dfloat);

  //compute the inverse diagonal
  for (dlong n=0;n<Nrows;n++) diagInv[n] = 1.0/diagA[n];
}

void parCSR::haloSetup(hlong *colIds) {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  hlong globalOffset = globalColStarts[rank];

  //collect the unique nonlocal column ids
  parallelId_t*  parIds = (parallelId_t*) malloc(offd->nnz*sizeof(parallelId_t));

  for (dlong n=0;n<offd->nnz;n++) {
    parIds[n].localId  = n;
    parIds[n].globalId = colIds[n];
  }

  //sort by global index
  qsort(parIds, offd->nnz, sizeof(parallelId_t), CompareGlobalId);

  //count unique nonlocal column ids
  dlong Noffdcols = 0; //number of unique columns
  if(offd->nnz) parIds[0].newId = Noffdcols;
  for (dlong n=1;n<offd->nnz;n++) {
    if (parIds[n].globalId != parIds[n-1].globalId)
      Noffdcols++;

    parIds[n].newId = Noffdcols;
  }
  if(offd->nnz) Noffdcols++;

  //record the global ids of the unique columns
  hlong *offdcols = (hlong *) malloc(Noffdcols*sizeof(hlong));
  Noffdcols = 0;
  if(offd->nnz) offdcols[Noffdcols++] = parIds[0].globalId;
  for (dlong n=1;n<offd->nnz;n++)
    if (parIds[n].globalId != parIds[n-1].globalId)
      offdcols[Noffdcols++] = parIds[n].globalId;

  //sort back to local order
  qsort(parIds, offd->nnz, sizeof(parallelId_t), CompareLocalId);

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
  for (dlong n=0;n<offd->nnz;n++)
    colIds[n] = NlocalCols + parIds[n].newId;

  free(parIds);
}

parCSR::~parCSR() {
  delete diag;
  delete offd;

  // free(diagA);
  // free(diagInv);

  if (o_diagA.size()) o_diagA.free();
  if (o_diagInv.size()) o_diagInv.free();

  // free(null);
  if (o_null.size()) o_null.free();

  // free(globalRowStarts);
  // free(globalColStarts);

  // free(colMap);

  if (halo)   halo->Free();
}

dfloat parCSR::rhoDinvA(){

  int rank, size;
  MPI_Comm_rank(comm, &rank);
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
  vectorRandomize(Nrows, Vx);

  dfloat norm_vo = vectorNorm(Nrows,Vx, comm);
  vectorScale(Nrows, 1.0/norm_vo, Vx);

  memcpy(V[0], Vx, Nrows*sizeof(dfloat));

  for(int j=0; j<k; j++){

    memcpy(Vx, V[j], Nrows*sizeof(dfloat));

    // v[j+1] = invD*(A*v[j])
    this->SpMV(1.0, Vx, 0., V[j+1]);
    vectorDotStar(Nrows, diagInv, V[j+1]);

    // modified Gram-Schmidth
    for(int i=0; i<=j; i++){
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = vectorInnerProd(Nrows, V[i], V[j+1],comm);

      // v[j+1] = v[j+1] - hij*v[i]
      vectorAdd(Nrows,-hij, V[i], 1.0, V[j+1]);

      H[i + j*k] = (double) hij;
    }

    if(j+1 < k){

      dfloat norm_vj = vectorNorm(Nrows,V[j+1],comm);

      H[j+1+ j*k] = (double) norm_vj;

      vectorScale(Nrows, 1./H[j+1 + j*k], V[j+1]);
    }
  }

  double *WR = (double *) calloc(k,sizeof(double));
  double *WI = (double *) calloc(k,sizeof(double));

  eig(k, H, WR, WI);

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

  occa::device& device = platform.device;

  diag->syncToDevice(device);
  offd->syncToDevice(device);

  if (Nrows) {
    o_diagA   = device.malloc(Nrows*sizeof(dfloat), diagA);
    o_diagInv = device.malloc(Nrows*sizeof(dfloat), diagInv);

    if(nullSpace)
      o_null = device.malloc(Nrows*sizeof(dfloat), null);
  }
}

} //namespace parAlmond
