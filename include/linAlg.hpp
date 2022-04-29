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

#ifndef LINALG_HPP
#define LINALG_HPP

#include "core.hpp"
#include "memory.hpp"

namespace libp {

class platform_t;

//launcher for basic linear algebra OCCA kernels
class linAlg_t {
 public:
  linAlg_t();
  linAlg_t(platform_t *_platform) { Setup(_platform); }

  void Setup(platform_t *_platform);

  //initialize list of kernels
  void InitKernels(std::vector<std::string> kernels);

  /*********************/
  /* vector operations */
  /*********************/

  // o_a[n] = alpha
  void set(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_a);

  // o_a[n] *= alpha
  void scale(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_a);

  // o_a[n] += alpha
  void add(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_a);

  // o_y[n] = beta*o_y[n] + alpha*o_x[n]
  void axpy(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_x,
                           const dfloat beta,  deviceMemory<dfloat> o_y);

  // o_z[n] = beta*o_y[n] + alpha*o_x[n]
  void zaxpy(const dlong N, const dfloat alpha, deviceMemory<dfloat> o_x,
                            const dfloat beta,  deviceMemory<dfloat> o_y,
                            deviceMemory<dfloat> o_z);

  // o_x[n] = alpha*o_a[n]*o_x[n]
  void amx(const dlong N, const dfloat alpha,
           deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x);

  // o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
  void amxpy(const dlong N, const dfloat alpha,
             deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
             const dfloat beta, deviceMemory<dfloat> o_y);

  // o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
  void zamxpy(const dlong N, const dfloat alpha,
              deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
              const dfloat beta, deviceMemory<dfloat> o_y, deviceMemory<dfloat> o_z);

  // o_x[n] = alpha*o_x[n]/o_a[n]
  void adx(const dlong N, const dfloat alpha,
           deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x);

  // o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  void adxpy(const dlong N, const dfloat alpha,
             deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
             const dfloat beta, deviceMemory<dfloat> o_y);

  // o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  void zadxpy(const dlong N, const dfloat alpha,
              deviceMemory<dfloat> o_a, deviceMemory<dfloat> o_x,
              const dfloat beta, deviceMemory<dfloat> o_y, deviceMemory<dfloat> o_z);

  // \min o_a
  dfloat min(const dlong N, deviceMemory<dfloat> o_a, comm_t comm);

  // \max o_a
  dfloat max(const dlong N, deviceMemory<dfloat> o_a, comm_t comm);

  // \sum o_a
  dfloat sum(const dlong N, deviceMemory<dfloat> o_a, comm_t comm);

  // ||o_a||_2
  dfloat norm2(const dlong N, deviceMemory<dfloat> o_a, comm_t comm);

  // o_x.o_y
  dfloat innerProd(const dlong N, deviceMemory<dfloat> o_x, deviceMemory<dfloat> o_y,
                    comm_t comm);

  // ||o_a||_w2
  dfloat weightedNorm2(const dlong N, deviceMemory<dfloat> o_w, deviceMemory<dfloat> o_a,
                       comm_t comm);

  // o_w.o_x.o_y
  dfloat weightedInnerProd(const dlong N, deviceMemory<dfloat> o_w, deviceMemory<dfloat> o_x,
                            deviceMemory<dfloat> o_y, comm_t comm);

  static void matrixRightSolve(const int NrowsA, const int NcolsA, const memory<double> A,
                               const int NrowsB, const int NcolsB, const memory<double> B,
                               memory<double> C);
  static void matrixRightSolve(const int NrowsA, const int NcolsA, const memory<float> A,
                               const int NrowsB, const int NcolsB, const memory<float> B,
                               memory<float> C);
  static void matrixUnderdeterminedRightSolveMinNorm(const int NrowsA, const int NcolsA,
                                                     const memory<double> A,
                                                     const memory<double> b,
                                                     memory<double> x);
  static void matrixUnderdeterminedRightSolveMinNorm(const int NrowsA, const int NcolsA,
                                                     const memory<float> A,
                                                     const memory<float> b,
                                                     memory<float> x);
  static void matrixUnderdeterminedRightSolveCPQR(const int NrowsA, const int NcolsA,
                                                  const memory<double> A,
                                                  const memory<double> b,
                                                  memory<double> x);
  static void matrixUnderdeterminedRightSolveCPQR(const int NrowsA, const int NcolsA,
                                                  const memory<float> A,
                                                  const memory<float> b,
                                                  memory<float> x);

  static void matrixEigenVectors(const int N, const memory<double> A,
                                 memory<double> VR, memory<double> WR, memory<double> WI);
  static void matrixEigenVectors(const int N, const memory<float> A,
                                 memory<float> VR, memory<float> WR, memory<float> WI);

  static void matrixEigenValues(const int N, const memory<double> A,
                                memory<double> WR, memory<double> WI);
  static void matrixEigenValues(const int N, const memory<float> A,
                                memory<float> WR, memory<float> WI);

  static void matrixInverse(const int N, memory<double> A);
  static void matrixInverse(const int N, memory<float> A);

  static double matrixConditionNumber(const int N, const memory<double> A);
  static float  matrixConditionNumber(const int N, const memory<float> A);

  static void matrixTranspose(const int M, const int N,
                             const memory<double> A, const int LDA,
                             memory<double> AT, const int LDAT);
  static void matrixTranspose(const int M, const int N,
                             const memory<float> A, const int LDA,
                             memory<float> AT, const int LDAT);

  static void matrixTranspose(const int M, const int N,
                              const memory<int> A, const int LDA,
                              memory<int> AT, const int LDAT);
  static void matrixTranspose(const int M, const int N,
                              const memory<long long int>  A, const int LDA,
                              memory<long long int> AT, const int LDAT);

 private:
  platform_t *platform;
  properties_t kernelInfo;

  static constexpr int blocksize = 256;

  //scratch space for reductions
  deviceMemory<dfloat> o_scratch;
  pinnedMemory<dfloat> h_scratch;

  kernel_t setKernel;
  kernel_t addKernel;
  kernel_t scaleKernel;
  kernel_t axpyKernel;
  kernel_t zaxpyKernel;
  kernel_t amxKernel;
  kernel_t amxpyKernel;
  kernel_t zamxpyKernel;
  kernel_t adxKernel;
  kernel_t adxpyKernel;
  kernel_t zadxpyKernel;
  kernel_t minKernel1;
  kernel_t minKernel2;
  kernel_t maxKernel1;
  kernel_t maxKernel2;
  kernel_t sumKernel1;
  kernel_t sumKernel2;
  kernel_t norm2Kernel1;
  kernel_t norm2Kernel2;
  kernel_t weightedNorm2Kernel1;
  kernel_t weightedNorm2Kernel2;
  kernel_t innerProdKernel1;
  kernel_t innerProdKernel2;
  kernel_t weightedInnerProdKernel1;
  kernel_t weightedInnerProdKernel2;
};

} //namespace libp

#endif
