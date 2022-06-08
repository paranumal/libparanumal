/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#ifndef PARALMOND_KERNELS_HPP
#define PARALMOND_KERNELS_HPP

namespace libp {

namespace parAlmond {

  void buildParAlmondKernels(platform_t& platform);

  void freeParAlmondKernels();

  //NC: Hard code these for now. Should be sufficient for GPU devices, but needs attention for CPU
  constexpr int blockSize = 256;
  constexpr int NonzerosPerBlock = 2048; //should be a multiple of blockSize for good unrolling

  extern kernel_t SpMVcsrKernel1;
  extern kernel_t SpMVcsrKernel2;
  extern kernel_t SpMVmcsrKernel;

  extern kernel_t SmoothJacobiCSRKernel;
  extern kernel_t SmoothJacobiMCSRKernel;

  extern kernel_t SmoothChebyshevStartKernel;
  extern kernel_t SmoothChebyshevCSRKernel;
  extern kernel_t SmoothChebyshevMCSRKernel;
  extern kernel_t SmoothChebyshevUpdateKernel;

  extern kernel_t vectorAddInnerProdKernel;
  extern kernel_t vectorAddWeightedInnerProdKernel;
  extern kernel_t kcycleCombinedOp1Kernel;
  extern kernel_t kcycleCombinedOp2Kernel;
  extern kernel_t kcycleWeightedCombinedOp1Kernel;
  extern kernel_t kcycleWeightedCombinedOp2Kernel;

  extern kernel_t dGEMVKernel;

} //namespace parAlmond

} // namespace libp

#endif
