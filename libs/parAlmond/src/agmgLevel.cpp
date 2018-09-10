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
#include "agmg.hpp"

namespace parAlmond {

agmgLevel::agmgLevel(parCSR *AA, KrylovType Ktype):
  multigridLevel(AA->Nrows, AA->Ncols, Ktype, AA->comm) {

  weighted = false;

  A = AA;
}

agmgLevel::agmgLevel(parCSR *AA, parCSR *PP, parCSR *RR, KrylovType Ktype):
  multigridLevel(AA->Nrows, AA->Ncols, Ktype, AA->comm) {

  //max
  Ncols = (AA->Ncols>PP->Ncols) ? AA->Ncols : PP->Ncols;

  weighted = false;

  A = AA;
  P = PP;
  R = RR;
}

agmgLevel::~agmgLevel() {

  delete   A; delete   P; delete   R;
  delete o_A; delete o_P; delete o_R;

}

void agmgLevel::Ax        (dfloat *x, dfloat *Ax){ A->SpMV(1.0, x, 0.0, Ax); }

void agmgLevel::coarsen   (dfloat *r, dfloat *Rr){ R->SpMV(1.0, r, 0.0, Rr); }

void agmgLevel::prolongate(dfloat *x, dfloat *Px){ P->SpMV(1.0, x, 1.0, Px); }

void agmgLevel::residual  (dfloat *rhs, dfloat *x, dfloat *res) { A->SpMV(-1.0, x, 1.0, rhs, res); }

void agmgLevel::Ax        (occa::memory o_x, occa::memory o_Ax){ o_A->SpMV(1.0, o_x, 0.0, o_Ax); }

void agmgLevel::coarsen   (occa::memory o_r, occa::memory o_Rr){ o_R->SpMV(1.0, o_r, 0.0, o_Rr); }

void agmgLevel::prolongate(occa::memory o_x, occa::memory o_Px){ o_P->SpMV(1.0, o_x, 1.0, o_Px); }

void agmgLevel::residual  (occa::memory o_rhs, occa::memory o_x, occa::memory o_res) { o_A->SpMV(-1.0, o_x, 1.0, o_rhs, o_res); }

void agmgLevel::smooth(dfloat *rhs, dfloat *x, bool x_is_zero){
  if(stype == JACOBI){
    this->smoothJacobi(rhs, x, x_is_zero);
  } else if(stype == DAMPED_JACOBI){
    this->smoothDampedJacobi(rhs, x, x_is_zero);
  } else if(stype == CHEBYSHEV){
    this->smoothChebyshev(rhs, x, x_is_zero);
  }
}

void agmgLevel::smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero){
  if(stype == JACOBI){
    this->smoothJacobi(o_rhs, o_x, x_is_zero);
  } else if(stype == DAMPED_JACOBI){
    this->smoothDampedJacobi(o_rhs, o_x, x_is_zero);
  } else if(stype == CHEBYSHEV){
    this->smoothChebyshev(o_rhs, o_x, x_is_zero);
  }
}

void agmgLevel::Report() {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

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
  dfloat avgNnz;
  MPI_Allreduce(&nnz, &maxNnz,   1, MPI_LONG_LONG_INT, MPI_MAX, comm);
  MPI_Allreduce(&nnz, &totalNnz, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);
  avgNnz = (dfloat) totalNnz/totalActive;

  if (nnz==0) nnz = maxNnz; //set this so it's ignored for the global min
  MPI_Allreduce(&nnz, &minNnz, 1, MPI_LONG_LONG_INT, MPI_MIN, comm);

  dfloat nnzPerRow = (Nrows==0) ? 0 : (dfloat) nnz/Nrows;
  dfloat minNnzPerRow=0, maxNnzPerRow=0, avgNnzPerRow=0;
  MPI_Allreduce(&nnzPerRow, &maxNnzPerRow, 1, MPI_DFLOAT, MPI_MAX, comm);
  MPI_Allreduce(&nnzPerRow, &avgNnzPerRow, 1, MPI_DFLOAT, MPI_SUM, comm);
  avgNnzPerRow /= totalActive;

  if (Nrows==0) nnzPerRow = maxNnzPerRow;
  MPI_Allreduce(&nnzPerRow, &minNnzPerRow, 1, MPI_DFLOAT, MPI_MIN, comm);

  if (rank==0){
    printf(     "|        %4d  |   %10.2f  |   %10.2f  |   %10.2f  |\n",
      totalActive, (dfloat)minNrows, (dfloat)minNnz, minNnzPerRow);
    printf("     |              |   %10.2f  |   %10.2f  |   %10.2f  |\n",
      (dfloat)maxNrows, (dfloat)maxNnz, maxNnzPerRow);
    printf("     |              |   %10.2f  |   %10.2f  |   %10.2f  |\n",
      avgNrows, avgNnz, avgNnzPerRow);
  }
}

} //namespace parAlmond