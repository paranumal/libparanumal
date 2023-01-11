/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
#include "platform.hpp"

namespace libp {

/*********************/
/* vector operations */
/*********************/

// o_a[n] = alpha
template <typename T>  
void linAlg_t<T>::set(const dlong N, const T alpha, deviceMemory<T> o_a) {
  setKernel(N, alpha, o_a);
}

// o_a[n] += alpha
template <typename T>  
void linAlg_t<T>::add(const dlong N, const T alpha, deviceMemory<T> o_a) {
  addKernel(N, alpha, o_a);
}

// o_a[n] *= alpha
template <typename T>  
void linAlg_t<T>::scale(const dlong N, const T alpha, deviceMemory<T> o_a)  {
  scaleKernel(N, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
template <typename T>    
void linAlg_t<T>::axpy(const dlong N, const T alpha, deviceMemory<T> o_x,
                    const T beta,  deviceMemory<T> o_y) {
  axpyKernel(N, alpha, o_x, beta, o_y);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
template <typename T>    
void linAlg_t<T>::zaxpy(const dlong N, const T alpha, deviceMemory<T> o_x,
                     const T beta, deviceMemory<T> o_y, deviceMemory<T> o_z) {
  zaxpyKernel(N, alpha, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_a[n]*o_x[n]
template <typename T>    
void linAlg_t<T>::amx(const dlong N, const T alpha,
                   deviceMemory<T> o_a, deviceMemory<T> o_x) {
  amxKernel(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
template <typename T>    
void linAlg_t<T>::amxpy(const dlong N, const T alpha,
                     deviceMemory<T> o_a, deviceMemory<T> o_x,
                     const T beta, deviceMemory<T> o_y) {
  amxpyKernel(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
template <typename T>    
void linAlg_t<T>::zamxpy(const dlong N, const T alpha,
                      deviceMemory<T> o_a, deviceMemory<T> o_x,
                      const T beta, deviceMemory<T> o_y, deviceMemory<T> o_z) {
  zamxpyKernel(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_x[n]/o_a[n]
template <typename T>    
void linAlg_t<T>::adx(const dlong N, const T alpha,
                   deviceMemory<T> o_a, deviceMemory<T> o_x) {
  adxKernel(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
template <typename T>    
void linAlg_t<T>::adxpy(const dlong N, const T alpha,
                     deviceMemory<T> o_a, deviceMemory<T> o_x,
                     const T beta, deviceMemory<T> o_y) {
  adxpyKernel(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  template <typename T>  
void linAlg_t<T>::zadxpy(const dlong N, const T alpha,
                      deviceMemory<T> o_a, deviceMemory<T> o_x,
                      const T beta, deviceMemory<T> o_y, deviceMemory<T> o_z) {
  zadxpyKernel(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// \min o_a
template <typename T>    
T linAlg_t<T>::min(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  //pinned scratch buffer
  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  minKernel(Nblock, N, o_a, o_scratch);

  T globalmin;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globalmin = h_scratch[0];
  } else {
    globalmin = std::numeric_limits<T>::max();
  }

  comm.Allreduce(globalmin, Comm::Min);

  return globalmin;
}

// \max o_a
template <typename T>    
T linAlg_t<T>::max(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  maxKernel(Nblock, N, o_a, o_scratch);

  T globalmax;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globalmax = h_scratch[0];
  } else {
    globalmax = -std::numeric_limits<T>::max();
  }

  comm.Allreduce(globalmax, Comm::Max);

  return globalmax;
}

// \sum o_a
template <typename T>    
T linAlg_t<T>::sum(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  sumKernel(Nblock, N, o_a, o_scratch);

  T globalsum;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globalsum = h_scratch[0];
  } else {
    globalsum = 0.0;
  }

  comm.Allreduce(globalsum, Comm::Sum);

  return globalsum;
}

// ||o_a||_2
template <typename T>    
T linAlg_t<T>::norm2(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  norm2Kernel(Nblock, N, o_a, o_scratch);

  T globalnorm;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globalnorm = h_scratch[0];
  } else {
    globalnorm = 0.0;
  }

  comm.Allreduce(globalnorm, Comm::Sum);

  return sqrt(globalnorm);
}

// o_x.o_y
template <typename T>    
T linAlg_t<T>::innerProd(const dlong N, deviceMemory<T> o_x, deviceMemory<T> o_y,
                           comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  innerProdKernel(Nblock, N, o_x, o_y, o_scratch);

  T globaldot;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globaldot = h_scratch[0];
  } else {
    globaldot = 0.0;
  }

  comm.Allreduce(globaldot, Comm::Sum);

  return globaldot;
}

// o_w.o_x.o_y
template <typename T>    
T linAlg_t<T>::weightedInnerProd(const dlong N, deviceMemory<T> o_w,
                                   deviceMemory<T> o_x, deviceMemory<T> o_y,
                                   comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  weightedInnerProdKernel(Nblock, N, o_w, o_x, o_y, o_scratch);

  T globaldot;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globaldot = h_scratch[0];
  } else {
    globaldot = 0.0;
  }

  comm.Allreduce(globaldot, Comm::Sum);

  return globaldot;
}

// ||o_a||_w2
template <typename T>    
T linAlg_t<T>::weightedNorm2(const dlong N, deviceMemory<T> o_w,
                               deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  weightedNorm2Kernel(Nblock, N, o_w, o_a, o_scratch);

  T globalnorm;
  if (Nblock>0) {
    h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
    platform->finish();
    globalnorm = h_scratch[0];
  } else {
    globalnorm = 0.0;
  }

  comm.Allreduce(globalnorm, Comm::Sum);

  return sqrt(globalnorm);
}

} //namespace libp
