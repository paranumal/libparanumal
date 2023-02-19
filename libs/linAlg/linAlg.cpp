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
template<>
void linAlg_t::set(const dlong N, const float alpha, deviceMemory<float> o_a) {
  setKernelFloat(N, alpha, o_a);
}

template<>
void linAlg_t::set(const dlong N, const double alpha, deviceMemory<double> o_a) {
  setKernelDouble(N, alpha, o_a);
}

// o_a[n] += alpha
template<>
void linAlg_t::add(const dlong N, const float alpha, deviceMemory<float> o_a) {
  addKernelFloat(N, alpha, o_a);
}
template<>
void linAlg_t::add(const dlong N, const double alpha, deviceMemory<double> o_a) {
  addKernelDouble(N, alpha, o_a);
}

// o_a[n] *= alpha
template<>
void linAlg_t::scale(const dlong N, const float alpha, deviceMemory<float> o_a)  {
  scaleKernelFloat(N, alpha, o_a);
}
template<>
void linAlg_t::scale(const dlong N, const double alpha, deviceMemory<double> o_a)  {
  scaleKernelDouble(N, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
template<>
void linAlg_t::axpy(const dlong N, const float alpha, deviceMemory<float> o_x,
                    const float beta,  deviceMemory<float> o_y) {
  axpyKernelFloat(N, alpha, o_x, beta, o_y);
}
template<>
void linAlg_t::axpy(const dlong N, const double alpha, deviceMemory<double> o_x,
                    const double beta,  deviceMemory<double> o_y) {
  axpyKernelDouble(N, alpha, o_x, beta, o_y);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
template<>
void linAlg_t::zaxpy(const dlong N, const float alpha, deviceMemory<float> o_x,
                     const float beta, deviceMemory<float> o_y, deviceMemory<float> o_z) {
  zaxpyKernelFloat(N, alpha, o_x, beta, o_y, o_z);
}
template<>
void linAlg_t::zaxpy(const dlong N, const double alpha, deviceMemory<double> o_x,
                     const double beta, deviceMemory<double> o_y, deviceMemory<double> o_z) {
  zaxpyKernelDouble(N, alpha, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_a[n]*o_x[n]
template<>
void linAlg_t::amx(const dlong N, const float alpha,
                   deviceMemory<float> o_a, deviceMemory<float> o_x) {
  amxKernelFloat(N, alpha, o_a, o_x);
}
template<>
void linAlg_t::amx(const dlong N, const double alpha,
                   deviceMemory<double> o_a, deviceMemory<double> o_x) {
  amxKernelDouble(N, alpha, o_a, o_x);
}


// o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
template<>
void linAlg_t::amxpy(const dlong N, const float alpha,
                     deviceMemory<float> o_a, deviceMemory<float> o_x,
                     const float beta, deviceMemory<float> o_y) {
  amxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y);
}
template<>
void linAlg_t::amxpy(const dlong N, const double alpha,
                     deviceMemory<double> o_a, deviceMemory<double> o_x,
                    const double beta, deviceMemory<double> o_y) {
  amxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
template<>
void linAlg_t::zamxpy(const dlong N, const float alpha,
                      deviceMemory<float> o_a, deviceMemory<float> o_x,
                      const float beta, deviceMemory<float> o_y, deviceMemory<float> o_z) {
  zamxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y, o_z);
}
template<>
void linAlg_t::zamxpy(const dlong N, const double alpha,
                      deviceMemory<double> o_a, deviceMemory<double> o_x,
                      const double beta, deviceMemory<double> o_y, deviceMemory<double> o_z) {
  zamxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_x[n]/o_a[n]
template<>
void linAlg_t::adx(const dlong N, const float alpha,
                   deviceMemory<float> o_a, deviceMemory<float> o_x) {
  adxKernelFloat(N, alpha, o_a, o_x);
}
template<>
void linAlg_t::adx(const dlong N, const double alpha,
                   deviceMemory<double> o_a, deviceMemory<double> o_x) {
  adxKernelDouble(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
template<>
void linAlg_t::adxpy(const dlong N, const float alpha,
                     deviceMemory<float> o_a, deviceMemory<float> o_x,
                     const float beta, deviceMemory<float> o_y) {
  adxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y);
}
template<>
void linAlg_t::adxpy(const dlong N, const double alpha,
                     deviceMemory<double> o_a, deviceMemory<double> o_x,
                     const double beta, deviceMemory<double> o_y) {
    adxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y);
}


// o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
template<>
void linAlg_t::zadxpy(const dlong N, const float alpha,
                      deviceMemory<float> o_a, deviceMemory<float> o_x,
                      const float beta, deviceMemory<float> o_y, deviceMemory<float> o_z) {
  zadxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y, o_z);
}
template<>
void linAlg_t::zadxpy(const dlong N, const double alpha,
                      deviceMemory<double> o_a, deviceMemory<double> o_x,
                     const double beta, deviceMemory<double> o_y, deviceMemory<double> o_z) {
  zadxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// \min o_a
template<typename T>
T linAlg_t::min(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  //pinned scratch buffer
  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    minKernelFloat(Nblock, N, o_a, o_scratch);
  } else {
    minKernelDouble(Nblock, N, o_a, o_scratch);
  }

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

template
float linAlg_t::min(const dlong N, deviceMemory<float> o_a, comm_t comm);

template
double linAlg_t::min(const dlong N, deviceMemory<double> o_a, comm_t comm);

// \max o_a
template<typename T>
T linAlg_t::max(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    maxKernelFloat(Nblock, N, o_a, o_scratch);
  } else {
    maxKernelDouble(Nblock, N, o_a, o_scratch);
  }

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

template
float linAlg_t::max(const dlong N, deviceMemory<float> o_a, comm_t comm);

template
double linAlg_t::max(const dlong N, deviceMemory<double> o_a, comm_t comm);


// \sum o_a
template<typename T>
T linAlg_t::sum(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    sumKernelFloat(Nblock, N, o_a, o_scratch);
  } else {
    sumKernelDouble(Nblock, N, o_a, o_scratch);
  }

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

template
float linAlg_t::sum(const dlong N, deviceMemory<float> o_a, comm_t comm);

template
double linAlg_t::sum(const dlong N, deviceMemory<double> o_a, comm_t comm);

// ||o_a||_2
template<typename T>
T linAlg_t::norm2(const dlong N, deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    norm2KernelFloat(Nblock, N, o_a, o_scratch);
  } else {
    norm2KernelDouble(Nblock, N, o_a, o_scratch);
  }

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

template
float linAlg_t::norm2(const dlong N, deviceMemory<float> o_a, comm_t comm);

template
double linAlg_t::norm2(const dlong N, deviceMemory<double> o_a, comm_t comm);

// o_x.o_y
template<typename T>
T linAlg_t::innerProd(const dlong N, deviceMemory<T> o_x, deviceMemory<T> o_y,
                         comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    innerProdKernelFloat(Nblock, N, o_x, o_y, o_scratch);
  } else {
    innerProdKernelDouble(Nblock, N, o_x, o_y, o_scratch);
  }

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

template
float linAlg_t::innerProd(const dlong N, deviceMemory<float> o_x, deviceMemory<float> o_y,
                         comm_t comm);

template
double linAlg_t::innerProd(const dlong N, deviceMemory<double> o_x, deviceMemory<double> o_y,
                         comm_t comm);

// o_w.o_x.o_y
template<typename T>
T linAlg_t::weightedInnerProd(const dlong N, deviceMemory<T> o_w,
                                  deviceMemory<T> o_x, deviceMemory<T> o_y,
                                  comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    weightedInnerProdKernelFloat(Nblock, N, o_w, o_x, o_y, o_scratch);
  } else {
    weightedInnerProdKernelDouble(Nblock, N, o_w, o_x, o_y, o_scratch);
  }

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

template
float linAlg_t::weightedInnerProd(const dlong N, deviceMemory<float> o_w,
                                  deviceMemory<float> o_x, deviceMemory<float> o_y,
                                  comm_t comm);

template
double linAlg_t::weightedInnerProd(const dlong N, deviceMemory<double> o_w,
                                  deviceMemory<double> o_x, deviceMemory<double> o_y,
                                  comm_t comm);

// ||o_a||_w2
template<typename T>
T linAlg_t::weightedNorm2(const dlong N, deviceMemory<T> o_w,
                              deviceMemory<T> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  pinnedMemory<T> h_scratch = platform->hostReserve<T>(1);
  deviceMemory<T> o_scratch = platform->reserve<T>(blocksize);

  if constexpr (std::is_same_v<T,float>) {
    weightedNorm2KernelFloat(Nblock, N, o_w, o_a, o_scratch);
  } else {
    weightedNorm2KernelDouble(Nblock, N, o_w, o_a, o_scratch);
  }

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

template
float linAlg_t::weightedNorm2(const dlong N, deviceMemory<float> o_w,
                              deviceMemory<float> o_a, comm_t comm);

template
double linAlg_t::weightedNorm2(const dlong N, deviceMemory<double> o_w,
                              deviceMemory<double> o_a, comm_t comm);

void linAlg_t::p2d(const dlong N, deviceMemory<pfloat> o_p, deviceMemory<dfloat> o_d) {
  p2dKernel(N, o_p, o_d);
}

void linAlg_t::d2p(const dlong N, deviceMemory<dfloat> o_d, deviceMemory<pfloat> o_p) {
  d2pKernel(N, o_d, o_p);
}


} //namespace libp
