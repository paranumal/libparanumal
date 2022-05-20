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
void linAlg_t::set(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_a) {
  setKernel(N, alpha, o_a);
}

// o_a[n] += alpha
void linAlg_t::add(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_a) {
  addKernel(N, alpha, o_a);
}

// o_a[n] *= alpha
void linAlg_t::scale(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_a)  {
  scaleKernel(N, alpha, o_a);
}

// o_y[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::axpy(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_x,
                    const dfloat beta,  deviceMemory<dfloat> o_y) {
  axpyKernel(N, alpha, o_x, beta, o_y);
}

// o_z[n] = beta*o_y[n] + alpha*o_x[n]
void linAlg_t::zaxpy(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_x,
                     const dfloat beta, deviceMemory<dfloat> o_y, deviceMemory<dfloat> o_z) {
  zaxpyKernel(N, alpha, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_a[n]*o_x[n]
void linAlg_t::amx(const dlong N, const dfloat alpha,
                   deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x) {
  amxKernel(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
void linAlg_t::amxpy(const dlong N, const dfloat alpha,
                     deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
                     const dfloat beta, deviceMemory<dfloat> o_y) {
  amxpyKernel(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
void linAlg_t::zamxpy(const dlong N, const dfloat alpha,
                      deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
                      const dfloat beta, deviceMemory<dfloat> o_y, deviceMemory<dfloat> o_z) {
  zamxpyKernel(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// o_x[n] = alpha*o_x[n]/o_a[n]
void linAlg_t::adx(const dlong N, const dfloat alpha,
                   deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x) {
  adxKernel(N, alpha, o_a, o_x);
}

// o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
void linAlg_t::adxpy(const dlong N, const dfloat alpha,
                     deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
                     const dfloat beta, deviceMemory<dfloat> o_y) {
  adxpyKernel(N, alpha, o_a, o_x, beta, o_y);
}

// o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
void linAlg_t::zadxpy(const dlong N, const dfloat alpha,
                      deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
                      const dfloat beta, deviceMemory<dfloat> o_y, deviceMemory<dfloat> o_z) {
  zadxpyKernel(N, alpha, o_a, o_x, beta, o_y, o_z);
}

// \min o_a
dfloat linAlg_t::min(const dlong N, deviceMemory<dfloat> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  minKernel1(Nblock, N, o_a, o_scratch);
  minKernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globalmin = h_scratch[0];
  comm.Allreduce(globalmin, Comm::Min);

  return globalmin;
}

// \max o_a
dfloat linAlg_t::max(const dlong N, deviceMemory<dfloat> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  maxKernel1(Nblock, N, o_a, o_scratch);
  maxKernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globalmax = h_scratch[0];
  comm.Allreduce(globalmax, Comm::Max);

  return globalmax;
}

// \sum o_a
dfloat linAlg_t::sum(const dlong N, deviceMemory<dfloat> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  sumKernel1(Nblock, N, o_a, o_scratch);
  sumKernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globalsum = h_scratch[0];
  comm.Allreduce(globalsum, Comm::Sum);

  return globalsum;
}

// ||o_a||_2
dfloat linAlg_t::norm2(const dlong N, deviceMemory<dfloat> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  norm2Kernel1(Nblock, N, o_a, o_scratch);
  norm2Kernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globalnorm = h_scratch[0];
  comm.Allreduce(globalnorm, Comm::Sum);

  return sqrt(globalnorm);
}

// o_x.o_y
dfloat linAlg_t::innerProd(const dlong N, deviceMemory<dfloat> o_x, deviceMemory<dfloat> o_y,
                           comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  innerProdKernel1(Nblock, N, o_x, o_y, o_scratch);
  innerProdKernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globaldot = h_scratch[0];
  comm.Allreduce(globaldot, Comm::Sum);

  return globaldot;
}

// o_w.o_x.o_y
dfloat linAlg_t::weightedInnerProd(const dlong N, deviceMemory<dfloat> o_w,
                                   deviceMemory<dfloat> o_x, deviceMemory<dfloat> o_y,
                                   comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  weightedInnerProdKernel1(Nblock, N, o_w, o_x, o_y, o_scratch);
  weightedInnerProdKernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globaldot = h_scratch[0];
  comm.Allreduce(globaldot, Comm::Sum);

  return globaldot;
}

// ||o_a||_w2
dfloat linAlg_t::weightedNorm2(const dlong N, deviceMemory<dfloat> o_w,
                               deviceMemory<dfloat> o_a, comm_t comm) {
  int Nblock = (N+blocksize-1)/blocksize;
  Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

  weightedNorm2Kernel1(Nblock, N, o_w, o_a, o_scratch);
  weightedNorm2Kernel2(Nblock, o_scratch);

  h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
  platform->finish();

  dfloat globalnorm = h_scratch[0];
  comm.Allreduce(globalnorm, Comm::Sum);

  return sqrt(globalnorm);
}

} //namespace libp
