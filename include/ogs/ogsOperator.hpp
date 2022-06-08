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

#ifndef OGS_OPERATOR_HPP
#define OGS_OPERATOR_HPP

#include "ogs.hpp"

namespace libp {

namespace ogs {

// The Z operator class is essentially a sparse CSR matrix,
// with no vals stored. By construction, the sparse
// matrix will have at most 1 non-zero per column.
class ogsOperator_t {
public:
  platform_t platform;

  dlong Ncols=0;
  dlong NrowsN=0;
  dlong NrowsT=0;
  dlong nnzN=0;
  dlong nnzT=0;

  memory<dlong> rowStartsN;
  memory<dlong> rowStartsT;
  memory<dlong> colIdsN;
  memory<dlong> colIdsT;

  deviceMemory<dlong> o_rowStartsN;
  deviceMemory<dlong> o_rowStartsT;
  deviceMemory<dlong> o_colIdsN;
  deviceMemory<dlong> o_colIdsT;

  dlong NrowBlocksN=0;
  dlong NrowBlocksT=0;
  memory<dlong> blockRowStartsN;
  memory<dlong> blockRowStartsT;
  deviceMemory<dlong> o_blockRowStartsN;
  deviceMemory<dlong> o_blockRowStartsT;

  Kind kind;

  ogsOperator_t()=default;
  ogsOperator_t(platform_t& _platform)
   : platform(_platform) {};

  void Free();

  void setupRowBlocks();

  //Apply Z operator
  template<template<typename> class U,
           template<typename> class V,
           typename T>
  void Gather(U<T> gv, const V<T> v,
              const int k, const Op op, const Transpose trans);

  template<typename T>
  void Gather(deviceMemory<T> gv, const deviceMemory<T> v,
              const int k, const Op op, const Transpose trans);

  //Apply Z^T transpose operator
  template<template<typename> class U,
           template<typename> class V,
           typename T>
  void Scatter(U<T> v, const V<T> gv,
               const int k, const Transpose trans);

  template<typename T>
  void Scatter(deviceMemory<T> v, const deviceMemory<T> gv,
               const int k, const Transpose trans);

  //Apply Z^T*Z operator
  template<template<typename> class U,
           typename T>
  void GatherScatter(U<T> v, const int k,
                     const Op op, const Transpose trans);

  template<typename T>
  void GatherScatter(deviceMemory<T> v, const int k,
                     const Op op, const Transpose trans);

private:
  template <template<typename> class U,
            template<typename> class V,
            template<typename> class Op,
            typename T>
  void Gather(U<T> gv, const V<T> v,
              const int K, const Transpose trans);
  template <template<typename> class U,
            template<typename> class Op,
            typename T>
  void GatherScatter(U<T> v, const int K,
                     const Transpose trans);

  //NC: Hard code these for now. Should be sufficient for GPU devices, but needs attention for CPU
  static constexpr int blockSize = 256;
  static constexpr int gatherNodesPerBlock = 512; //should be a multiple of blockSize for good unrolling

  //4 types - Float, Double, Int32, Int64
  //4 ops - Add, Mul, Max, Min
  static kernel_t gatherScatterKernel[4][4];
  static kernel_t gatherKernel[4][4];
  static kernel_t scatterKernel[4];

  friend void InitializeKernels(platform_t& platform, const Type type, const Op op);
};

template <template<typename> class U,
          template<typename> class V,
          typename T>
void extract(const dlong N,
             const int K,
             const memory<dlong> ids,
             const U<T> q,
             V<T> gatherq);

} //namespace ogs

} //namespace libp

#endif
