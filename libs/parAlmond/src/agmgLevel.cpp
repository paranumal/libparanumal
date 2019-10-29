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

agmgLevel::agmgLevel(parCSR *A_, KrylovType ktype_):
  multigridLevel(A_->Nrows, A_->Ncols, ktype_, A_->comm) {

  weighted = false;
  gatherLevel = false;
  ogs=NULL;

  A = A_;
  P = NULL;
  R = NULL;
  o_A = NULL;
  o_P = NULL;
  o_R = NULL;
}

agmgLevel::agmgLevel(parCSR *A_, parCSR *P_, parCSR *R_, KrylovType ktype_):
  multigridLevel(A_->Nrows, A_->Ncols, ktype_, A_->comm) {

  //max
  Ncols = (A_->Ncols>P_->Ncols) ? A_->Ncols : P_->Ncols;

  weighted = false;
  gatherLevel = false;
  ogs=NULL;

  A = A_;
  P = P_;
  R = R_;
  o_A = NULL;
  o_P = NULL;
  o_R = NULL;
}

agmgLevel::~agmgLevel() {

  if (  A) delete   A;
  if (  P) delete   P;
  if (  R) delete   R;
  if (o_A) delete o_A;
  if (o_P) delete o_P;
  if (o_R) delete o_R;
}

void agmgLevel::Ax        (dfloat *X, dfloat *Ax){ A->SpMV(1.0, X, 0.0, Ax); }

void agmgLevel::coarsen   (dfloat *r, dfloat *Rr){
  if (gatherLevel) {
    ogs->Gather(Gx, r, ogs_dfloat, ogs_add, ogs_notrans);
    R->SpMV(1.0, Gx, 0.0, Rr);
  } else {
    R->SpMV(1.0, r, 0.0, Rr);
  }
}

void agmgLevel::prolongate(dfloat *X, dfloat *Px){
  if (gatherLevel) {
    P->SpMV(1.0, X, 0.0, Gx);
    ogs->Scatter(Sx, Gx, ogs_dfloat, ogs_add, ogs_notrans);
    vectorAdd(P->Nrows, 1.0, Sx, 1.0, Px);
  } else {
    P->SpMV(1.0, X, 1.0, Px);
  }
}

void agmgLevel::residual  (dfloat *RHS, dfloat *X, dfloat *RES) { A->SpMV(-1.0, X, 1.0, RHS, RES); }

void agmgLevel::Ax        (occa::memory& o_X, occa::memory& o_Ax){ o_A->SpMV(1.0, o_X, 0.0, o_Ax); }

void agmgLevel::coarsen   (occa::memory& o_r, occa::memory& o_Rr){
  if (gatherLevel) {
    ogs->Gather(o_Gx, o_r, ogs_dfloat, ogs_add, ogs_notrans);
    // vectorDotStar(ogs->Ngather, o_gatherWeight, o_Gx);
    o_R->SpMV(1.0, o_Gx, 0.0, o_Rr);
  } else {
    o_R->SpMV(1.0, o_r, 0.0, o_Rr);
  }
}

void agmgLevel::prolongate(occa::memory& o_X, occa::memory& o_Px){
  if (gatherLevel) {
    o_P->SpMV(1.0, o_X, 0.0, o_Gx);
    ogs->Scatter(o_Sx, o_Gx, ogs_dfloat, ogs_add, ogs_notrans);
    vectorAdd(ogs->N, 1.0, o_Sx, 1.0, o_Px);
  } else {
    o_P->SpMV(1.0, o_X, 1.0, o_Px);
  }
}

void agmgLevel::residual  (occa::memory& o_RHS, occa::memory& o_X, occa::memory& o_RES) { o_A->SpMV(-1.0, o_X, 1.0, o_RHS, o_RES); }

void agmgLevel::smooth(dfloat *RHS, dfloat *X, bool x_is_zero){
  if(stype == JACOBI){
    this->smoothJacobi(RHS, X, x_is_zero);
  } else if(stype == DAMPED_JACOBI){
    this->smoothDampedJacobi(RHS, X, x_is_zero);
  } else if(stype == CHEBYSHEV){
    this->smoothChebyshev(RHS, X, x_is_zero);
  }
}

void agmgLevel::smooth(occa::memory& o_RHS, occa::memory& o_X, bool x_is_zero){
  if(stype == JACOBI){
    this->smoothJacobi(o_RHS, o_X, x_is_zero);
  } else if(stype == DAMPED_JACOBI){
    this->smoothDampedJacobi(o_RHS, o_X, x_is_zero);
  } else if(stype == CHEBYSHEV){
    this->smoothChebyshev(o_RHS, o_X, x_is_zero);
  }
}

void agmgLevel::Report() {

  //This setup can be called by many subcommunicators, so only
  // print on the global root.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  hlong hNrows = (hlong) Nrows;

  int active = (Nrows>0) ? 1:0;
  int totalActive=0;
  MPI_Allreduce(&active, &totalActive, 1, MPI_INT, MPI_SUM, comm);

  dlong minNrows=0, maxNrows=0;
  hlong totalNrows=0;
  dfloat avgNrows;
  MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, comm);
  avgNrows = (dfloat) totalNrows/totalActive;

  if (Nrows==0) Nrows=maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, comm);


  long long int nnz;
  nnz = A->diag->nnz+A->offd->nnz;

  long long int minNnz=0, maxNnz=0, totalNnz=0;
  MPI_Allreduce(&nnz, &maxNnz,   1, MPI_LONG_LONG_INT, MPI_MAX, comm);
  MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

  if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
  MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, comm);

  dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
  dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
  MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, comm);
  MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, comm);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) nnzPerRow = maxNnzPerRow;
  MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, comm);

  char smootherString[BUFSIZ];
  if (stype==DAMPED_JACOBI)
    strcpy(smootherString, "Damped Jacobi   ");
  else if (stype==CHEBYSHEV)
    strcpy(smootherString, "Chebyshev       ");

  if (rank==0){
    printf(     "|  parAlmond |  %12d  | %13d   |   %s|\n", minNrows, (int)minNnzPerRow, smootherString);
    printf("     |            |  %12d  | %13d   |                   |\n", maxNrows, (int)maxNnzPerRow);
    printf("     |            |  %12d  | %13d   |                   |\n", (int)avgNrows, (int)avgNnzPerRow);
  }
}

} //namespace parAlmond