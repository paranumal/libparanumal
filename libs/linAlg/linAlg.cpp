/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "linAlg.hpp"

/*********************/
/* vector operations */
/*********************/

// o_a[n] = alpha
template <class lafloat>
void linAlg_t<lafloat>::set(const dlong N, const lafloat alpha, occa::memory& o_a) {
  setKernel(N, alpha, o_a);
}

// o_a[n] += alpha
template <class lafloat>
void linAlg_t<lafloat>::add(const dlong N, const lafloat alpha, occa::memory& o_a) {
  addKernel(N, alpha, o_a);
}

// o_a[n] *= alpha
template <class lafloat>
void linAlg_t<lafloat>::scale(const dlong N, const lafloat alpha, occa::memory& o_a)  {
  scaleKernel(N, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
template <class lafloat>
void linAlg_t<lafloat>::axpy(const dlong N, const lafloat alpha, occa::memory& o_x,
                    const lafloat beta,  occa::memory& o_y) {
  axpyKernel(N, alpha, o_x, beta, o_y);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
template <class lafloat>
void linAlg_t<lafloat>::zaxpy(const dlong N, const lafloat alpha, occa::memory& o_x,
                     const lafloat beta, occa::memory& o_y, occa::memory& o_z) {
  zaxpyKernel(N, alpha, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_a[n]*o_x[n]
template <class lafloat>
void linAlg_t<lafloat>::amx(const dlong N, const lafloat alpha,
                   occa::memory& o_a, occa::memory& o_x) {
  amxKernel(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
template <class lafloat>
void linAlg_t<lafloat>::amxpy(const dlong N, const lafloat alpha,
                     occa::memory& o_a, occa::memory& o_x,
                     const lafloat beta, occa::memory& o_y) {
  amxpyKernel(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
template <class lafloat>
void linAlg_t<lafloat>::zamxpy(const dlong N, const lafloat alpha,
                      occa::memory& o_a, occa::memory& o_x,
                      const lafloat beta, occa::memory& o_y, occa::memory& o_z) {
  zamxpyKernel(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_x[n]/o_a[n]
template <class lafloat>
void linAlg_t<lafloat>::adx(const dlong N, const lafloat alpha,
                   occa::memory& o_a, occa::memory& o_x) {
  adxKernel(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
template <class lafloat>
void linAlg_t<lafloat>::adxpy(const dlong N, const lafloat alpha,
                     occa::memory& o_a, occa::memory& o_x,
                     const lafloat beta, occa::memory& o_y) {
  adxpyKernel(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
template <class lafloat>
void linAlg_t<lafloat>::zadxpy(const dlong N, const lafloat alpha,
                      occa::memory& o_a, occa::memory& o_x,
                      const lafloat beta, occa::memory& o_y, occa::memory& o_z) {
  zadxpyKernel(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// \min o_a
template <class lafloat>
lafloat linAlg_t<lafloat>::min(const dlong N, occa::memory& o_a, MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  minKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat min = std::numeric_limits<dfloat>::max();
  for(dlong n=0;n<Nblock;++n){
    min = (scratch[n]<min) ? scratch[n] : min;
  }

  dfloat globalmin = 0;
  MPI_Allreduce(&min, &globalmin, 1, MPI_DFLOAT, MPI_MIN, comm);

  return globalmin;
}

// \max o_a
template <class lafloat>
lafloat linAlg_t<lafloat>::max(const dlong N, occa::memory& o_a, MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  maxKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat max = -std::numeric_limits<dfloat>::max();
  for(dlong n=0;n<Nblock;++n){
    max = (scratch[n]>max) ? scratch[n] : max;
  }

  dfloat globalmax = 0;
  MPI_Allreduce(&max, &globalmax, 1, MPI_DFLOAT, MPI_MAX, comm);

  return globalmax;
}

// \sum o_a
template <class lafloat>
lafloat linAlg_t<lafloat>::sum(const dlong N, occa::memory& o_a, MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  sumKernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat sum = 0;
  for(dlong n=0;n<Nblock;++n){
    sum += scratch[n];
  }

  dfloat globalsum = 0;
  MPI_Allreduce(&sum, &globalsum, 1, MPI_DFLOAT, MPI_SUM, comm);

  return globalsum;
}

// ||o_a||_2
template <class lafloat>
lafloat linAlg_t<lafloat>::norm2(const dlong N, occa::memory& o_a, MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  norm2Kernel(Nblock, N, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat norm = 0;
  for(dlong n=0;n<Nblock;++n){
    norm += scratch[n];
  }

  dfloat globalnorm = 0;
  MPI_Allreduce(&norm, &globalnorm, 1, MPI_DFLOAT, MPI_SUM, comm);

  return sqrt(globalnorm);
}

// o_x.o_y
template <class lafloat>
lafloat linAlg_t<lafloat>::innerProd(const dlong N, occa::memory& o_x, occa::memory& o_y,
                           MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  innerProdKernel(Nblock, N, o_x, o_y, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat dot = 0;
  for(dlong n=0;n<Nblock;++n){
    dot += scratch[n];
  }

  dfloat globaldot = 0;
  MPI_Allreduce(&dot, &globaldot, 1, MPI_DFLOAT, MPI_SUM, comm);

  return globaldot;
}

// o_w.o_x.o_y
template <class lafloat>
lafloat linAlg_t<lafloat>::weightedInnerProd(const dlong N, occa::memory& o_w,
                                   occa::memory& o_x, occa::memory& o_y,
                                   MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  weightedInnerProdKernel(Nblock, N, o_w, o_x, o_y, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat dot = 0;
  for(dlong n=0;n<Nblock;++n){
    dot += scratch[n];
  }

  dfloat globaldot = 0;
  MPI_Allreduce(&dot, &globaldot, 1, MPI_DFLOAT, MPI_SUM, comm);

  return globaldot;
}

// ||o_a||_w2
template <class lafloat>
lafloat linAlg_t<lafloat>::weightedNorm2(const dlong N, occa::memory& o_w,
                               occa::memory& o_a, MPI_Comm comm) {
  //TODO, maybe complete reduction on device with second kernel?
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  weightedNorm2Kernel(Nblock, N, o_w, o_a, o_scratch);

  o_scratch.copyTo(scratch, Nblock*sizeof(lafloat));

  dfloat norm = 0;
  for(dlong n=0;n<Nblock;++n){
    norm += scratch[n];
  }

  dfloat globalnorm = 0;
  MPI_Allreduce(&norm, &globalnorm, 1, MPI_DFLOAT, MPI_SUM, comm);

  return sqrt(globalnorm);
}

template class linAlg_t<float>;
template class linAlg_t<double>;
