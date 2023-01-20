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
  template<typename T>
  void set(const dlong N, const T alpha, deviceMemory<T> o_a);

  // o_a[n] *= alpha
  template<typename T>
  void scale(const dlong N, const T alpha, deviceMemory<T> o_a);

  // o_a[n] += alpha
  template<typename T>
  void add(const dlong N, const T alpha, deviceMemory<T> o_a);

  // o_y[n] = beta*o_y[n] + alpha*o_x[n]
  template<typename T>
  void axpy(const dlong N, const T alpha, deviceMemory<T> o_x,
                           const T beta,  deviceMemory<T> o_y);

  // o_z[n] = beta*o_y[n] + alpha*o_x[n]
  template<typename T>
  void zaxpy(const dlong N, const T alpha, deviceMemory<T> o_x,
                            const T beta,  deviceMemory<T> o_y,
                            deviceMemory<T> o_z);

  // o_x[n] = alpha*o_a[n]*o_x[n]
  template<typename T>
  void amx(const dlong N, const T alpha,
           deviceMemory<T> o_a, deviceMemory<T> o_x);

  // o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
  template<typename T>
  void amxpy(const dlong N, const T alpha,
             deviceMemory<T> o_a, deviceMemory<T> o_x,
             const T beta, deviceMemory<T> o_y);

  // o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
  template<typename T>
  void zamxpy(const dlong N, const T alpha,
              deviceMemory<T> o_a, deviceMemory<T> o_x,
              const T beta, deviceMemory<T> o_y, deviceMemory<T> o_z);

  // o_x[n] = alpha*o_x[n]/o_a[n]
  template<typename T>
  void adx(const dlong N, const T alpha,
           deviceMemory<T> o_a, deviceMemory<T> o_x);

  // o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  template<typename T>
  void adxpy(const dlong N, const T alpha,
             deviceMemory<T> o_a, deviceMemory<T> o_x,
             const T beta, deviceMemory<T> o_y);

  // o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
  template<typename T>
  void zadxpy(const dlong N, const T alpha,
              deviceMemory<T> o_a, deviceMemory<T> o_x,
              const T beta, deviceMemory<T> o_y, deviceMemory<T> o_z);

  // \min o_a
  template<typename T>
  T min(const dlong N, deviceMemory<T> o_a, comm_t comm);

  // \max o_a
  template<typename T>
  T max(const dlong N, deviceMemory<T> o_a, comm_t comm);

  // \sum o_a
  template<typename T>
  T sum(const dlong N, deviceMemory<T> o_a, comm_t comm);

  // ||o_a||_2
  template<typename T>
  T norm2(const dlong N, deviceMemory<T> o_a, comm_t comm);

  // o_x.o_y
  template<typename T>
  T innerProd(const dlong N, deviceMemory<T> o_x, deviceMemory<T> o_y,
                    comm_t comm);

  // ||o_a||_w2
  template<typename T>
  T weightedNorm2(const dlong N, deviceMemory<T> o_w, deviceMemory<T> o_a,
                       comm_t comm);

  // o_w.o_x.o_y
  template<typename T>
  T weightedInnerProd(const dlong N, deviceMemory<T> o_w, deviceMemory<T> o_x,
                            deviceMemory<T> o_y, comm_t comm);

  // p=d
  void p2d(const dlong N, deviceMemory<pfloat> o_p, deviceMemory<dfloat> o_d);
  void d2p(const dlong N, deviceMemory<dfloat> o_d, deviceMemory<pfloat> o_p);

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
  properties_t kernelInfoFloat;
  properties_t kernelInfoDouble;

  static constexpr int blocksize = 256;

  kernel_t setKernelFloat;
  kernel_t addKernelFloat;
  kernel_t scaleKernelFloat;
  kernel_t axpyKernelFloat;
  kernel_t zaxpyKernelFloat;
  kernel_t amxKernelFloat;
  kernel_t amxpyKernelFloat;
  kernel_t zamxpyKernelFloat;
  kernel_t adxKernelFloat;
  kernel_t adxpyKernelFloat;
  kernel_t zadxpyKernelFloat;
  kernel_t minKernelFloat;
  kernel_t maxKernelFloat;
  kernel_t sumKernelFloat;
  kernel_t norm2KernelFloat;
  kernel_t weightedNorm2KernelFloat;
  kernel_t innerProdKernelFloat;
  kernel_t weightedInnerProdKernelFloat;

  kernel_t setKernelDouble;
  kernel_t addKernelDouble;
  kernel_t scaleKernelDouble;
  kernel_t axpyKernelDouble;
  kernel_t zaxpyKernelDouble;
  kernel_t amxKernelDouble;
  kernel_t amxpyKernelDouble;
  kernel_t zamxpyKernelDouble;
  kernel_t adxKernelDouble;
  kernel_t adxpyKernelDouble;
  kernel_t zadxpyKernelDouble;
  kernel_t minKernelDouble;
  kernel_t maxKernelDouble;
  kernel_t sumKernelDouble;
  kernel_t norm2KernelDouble;
  kernel_t weightedNorm2KernelDouble;
  kernel_t innerProdKernelDouble;
  kernel_t weightedInnerProdKernelDouble;

  kernel_t p2dKernel;
  kernel_t d2pKernel;

};


} //namespace libp

#endif
