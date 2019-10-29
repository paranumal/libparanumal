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

#ifndef LINALG_HPP
#define LINALG_HPP

#include <occa.hpp>
#include "types.h"
#include "utils.hpp"
#include "settings.hpp"

//launcher for basic linear algebra OCCA kernels
class linAlg_t {
public:
  occa::device& device;
  settings_t& settings;
  occa::properties& props;

  int blocksize;

  //scratch space for reductions
  dfloat *scratch;
  occa::memory h_scratch;
  occa::memory o_scratch;

  linAlg_t(occa::device& device_,
           settings_t& settings_, occa::properties& props_);

  //named constructor
  static linAlg_t& Setup(occa::device& device_,
           settings_t& settings_, occa::properties& props_);

  //initialize list of kernels
  void InitKernels(vector<string> kernels, MPI_Comm& comm);

  ~linAlg_t();

  /*********************/
  /* vector operations */
  /*********************/

  // o_a[n] = alpha
  void set(const dlong N, const dfloat alpha, occa::memory& o_a);

  // o_a[n] *= alpha
  void scale(const dlong N, const dfloat alpha, occa::memory& o_a);

  // o_a[n] += alpha
  void add(const dlong N, const dfloat alpha, occa::memory& o_a);

  // o_y[n] = beta*o_y[n] + alpha*o_x[n]
  void axpy(const dlong N, const dfloat alpha, occa::memory& o_x,
                           const dfloat beta,  occa::memory& o_y);

  // o_z[n] = beta*o_y[n] + alpha*o_x[n]
  void zaxpy(const dlong N, const dfloat alpha, occa::memory& o_x,
                            const dfloat beta,  occa::memory& o_y,
                            occa::memory& o_z);

  // o_x[n] = alpha*o_a[n]*o_x[n]
  void amx(const dlong N, const dfloat alpha,
           occa::memory& o_a, occa::memory& o_x);

  // o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
  void amxpy(const dlong N, const dfloat alpha,
             occa::memory& o_a, occa::memory& o_x,
             const dfloat beta, occa::memory& o_y);

  // o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
  void zamxpy(const dlong N, const dfloat alpha,
              occa::memory& o_a, occa::memory& o_x,
              const dfloat beta, occa::memory& o_y, occa::memory& o_z);

  // o_x[n] = alpha*o_x[n]/o_a[n]
  void adx(const dlong N, const dfloat alpha,
           occa::memory& o_a, occa::memory& o_x);

  // o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  void adxpy(const dlong N, const dfloat alpha,
             occa::memory& o_a, occa::memory& o_x,
             const dfloat beta, occa::memory& o_y);

  // o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  void zadxpy(const dlong N, const dfloat alpha,
              occa::memory& o_a, occa::memory& o_x,
              const dfloat beta, occa::memory& o_y, occa::memory& o_z);

  // \sum o_a
  dfloat sum(const dlong N, occa::memory& o_a, MPI_Comm comm);

  // ||o_a||_2
  dfloat norm2(const dlong N, occa::memory& o_a, MPI_Comm comm);

  // o_x.o_y
  dfloat innerProd(const dlong N, occa::memory& o_x, occa::memory& o_y,
                    MPI_Comm comm);

  // ||o_a||_w2
  dfloat weightedNorm2(const dlong N, occa::memory& o_w, occa::memory& o_a,
                       MPI_Comm comm);

  // o_w.o_x.o_y
  dfloat weightedInnerProd(const dlong N, occa::memory& o_w, occa::memory& o_x,
                            occa::memory& o_y, MPI_Comm comm);

  occa::kernel setKernel;
  occa::kernel addKernel;
  occa::kernel scaleKernel;
  occa::kernel axpyKernel;
  occa::kernel zaxpyKernel;
  occa::kernel amxKernel;
  occa::kernel amxpyKernel;
  occa::kernel zamxpyKernel;
  occa::kernel adxKernel;
  occa::kernel adxpyKernel;
  occa::kernel zadxpyKernel;
  occa::kernel sumKernel;
  occa::kernel norm2Kernel;
  occa::kernel weightedNorm2Kernel;
  occa::kernel innerProdKernel;
  occa::kernel weightedInnerProdKernel;
};

#endif
