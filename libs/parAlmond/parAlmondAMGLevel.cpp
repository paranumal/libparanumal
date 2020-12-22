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
#include "parAlmond/parAlmondMultigrid.hpp"
#include "parAlmond/parAlmondAMGLevel.hpp"

namespace parAlmond {

amgLevel::amgLevel(parCSR *_A, settings_t& _settings):
  multigridLevel(_A->Nrows, _A->Ncols, _A->platform, _settings),
  A(_A) {

  //determine smoother
  if (settings.compareSetting("PARALMOND SMOOTHER", "CHEBYSHEV")) {
    stype = CHEBYSHEV;
    settings.getSetting("PARALMOND CHEBYSHEV DEGREE", ChebyshevIterations);
  } else { //default to DAMPED_JACOBI
    stype = DAMPED_JACOBI;
  }
}

amgLevel::~amgLevel() {
  if (  A) delete   A;
  if (  P) delete   P;
  if (  R) delete   R;
}

void amgLevel::Operator(occa::memory& o_X, occa::memory& o_Ax){
  A->SpMV(1.0, o_X, 0.0, o_Ax);
}

void amgLevel::coarsen   (occa::memory& o_r, occa::memory& o_Rr){
  if (gatherLevel) {
    ogs->Gather(o_Gx, o_r, ogs_dfloat, ogs_add, ogs_notrans);
    // vectorDotStar(ogs->Ngather, o_gatherWeight, o_Gx);
    R->SpMV(1.0, o_Gx, 0.0, o_Rr);
  } else {
    R->SpMV(1.0, o_r, 0.0, o_Rr);
  }
}

void amgLevel::prolongate(occa::memory& o_X, occa::memory& o_Px){
  if (gatherLevel) {
    P->SpMV(1.0, o_X, 0.0, o_Gx);
    ogs->Scatter(o_Sx, o_Gx, ogs_dfloat, ogs_add, ogs_notrans);
    platform.linAlg.axpy(ogs->N, 1.0, o_Sx, 1.0, o_Px);
  } else {
    P->SpMV(1.0, o_X, 1.0, o_Px);
  }
}

void amgLevel::residual  (occa::memory& o_RHS, occa::memory& o_X,
                          occa::memory& o_RES) {
  A->SpMV(-1.0, o_X, 1.0, o_RHS, o_RES);
}

void amgLevel::smooth(occa::memory& o_RHS, occa::memory& o_X, bool x_is_zero){
  if(stype == DAMPED_JACOBI){
    A->smoothDampedJacobi(o_RHS, o_X, lambda,
                          x_is_zero, o_scratch);
  } else if(stype == CHEBYSHEV){
    A->smoothChebyshev(o_RHS, o_X, lambda0, lambda1,
                       x_is_zero, o_scratch,
                       ChebyshevIterations);
  }
}

void amgLevel::setupSmoother(){

  if (stype == DAMPED_JACOBI) {
    lambda = (4./3.)/A->rho;
  } else if (stype == CHEBYSHEV) {
    lambda1 = A->rho;
    lambda0 = A->rho/10.;
  }
}

void amgLevel::syncToDevice(){
  A->syncToDevice();
  if (P) P->syncToDevice();
  if (R) R->syncToDevice();
}

void amgLevel::Report() {

  //This setup can be called by many subcommunicators, so only
  // print on the global root.
  int rank;
  MPI_Comm_rank(A->comm, &rank);

  hlong hNrows = (hlong) Nrows;

  int active = (Nrows>0) ? 1:0;
  int totalActive=0;
  MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, A->comm);

  dlong minNrows=0, maxNrows=0;
  hlong totalNrows=0;
  dfloat avgNrows;
  MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, A->comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, A->comm);
  avgNrows = (dfloat) totalNrows/totalActive;

  if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, A->comm);


  long long int nnz;
  nnz = A->diag.nnz+A->offd.nnz;

  long long int minNnz=0, maxNnz=0, totalNnz=0;
  MPI_Allreduce(&nnz, &maxNnz,   1, MPI_LONG_LONG_INT, MPI_MAX, A->comm);
  MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, A->comm);

  if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
  MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, A->comm);

  dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
  dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
  MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, A->comm);
  MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, A->comm);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) nnzPerRow = maxNnzPerRow;
  MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, A->comm);

  char smootherString[BUFSIZ];
  if (stype==DAMPED_JACOBI)
    strcpy(smootherString, "Damped Jacobi   ");
  else if (stype==CHEBYSHEV)
    strcpy(smootherString, "Chebyshev       ");

  if (rank==0){
    printf(      "|  parAlmond |  %12lld  |  %12d  | %13d   |   %s|\n", (long long int) totalNrows, minNrows, (int)minNnzPerRow, smootherString);
    printf("      |            |                |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("      |            |                |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

} //namespace parAlmond