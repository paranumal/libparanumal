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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondKernels.hpp"

namespace libp {

namespace parAlmond {

kernel_t SpMVcsrKernel1;
kernel_t SpMVcsrKernel2;
kernel_t SpMVmcsrKernel;

kernel_t SmoothJacobiCSRKernel;
kernel_t SmoothJacobiMCSRKernel;

kernel_t SmoothChebyshevStartKernel;
kernel_t SmoothChebyshevCSRKernel;
kernel_t SmoothChebyshevMCSRKernel;
kernel_t SmoothChebyshevUpdateKernel;

kernel_t kcycleCombinedOp1Kernel;
kernel_t kcycleCombinedOp2Kernel;
kernel_t vectorAddInnerProdKernel;

kernel_t dGEMVKernel;

void buildParAlmondKernels(platform_t& platform){

  if (SpMVcsrKernel1.isInitialized()==false) {
    //seed rng
    int rank=platform.rank();
    double seed = (double) rank;
    srand48(seed);

    //build kernels
    properties_t kernelInfo = platform.props();

    kernelInfo["defines/" "p_BLOCKSIZE"]= blockSize;
    kernelInfo["defines/" "p_NonzerosPerBlock"]= NonzerosPerBlock;

    if (rank==0) {printf("Compiling parALMOND Kernels...");fflush(stdout);}

    SpMVcsrKernel1  = platform.buildKernel(PARALMOND_DIR"/okl/SpMVcsr.okl",  "SpMVcsr1",  kernelInfo);
    SpMVcsrKernel2  = platform.buildKernel(PARALMOND_DIR"/okl/SpMVcsr.okl",  "SpMVcsr2",  kernelInfo);
    SpMVmcsrKernel  = platform.buildKernel(PARALMOND_DIR"/okl/SpMVmcsr.okl", "SpMVmcsr1", kernelInfo);

    SmoothJacobiCSRKernel  = platform.buildKernel(PARALMOND_DIR"/okl/SmoothJacobi.okl", "SmoothJacobiCSR", kernelInfo);
    SmoothJacobiMCSRKernel = platform.buildKernel(PARALMOND_DIR"/okl/SmoothJacobi.okl", "SmoothJacobiMCSR", kernelInfo);

    SmoothChebyshevStartKernel = platform.buildKernel(PARALMOND_DIR"/okl/SmoothChebyshev.okl", "SmoothChebyshevStart", kernelInfo);
    SmoothChebyshevCSRKernel  = platform.buildKernel(PARALMOND_DIR"/okl/SmoothChebyshev.okl", "SmoothChebyshevCSR", kernelInfo);
    SmoothChebyshevMCSRKernel = platform.buildKernel(PARALMOND_DIR"/okl/SmoothChebyshev.okl", "SmoothChebyshevMCSR", kernelInfo);
    SmoothChebyshevUpdateKernel = platform.buildKernel(PARALMOND_DIR"/okl/SmoothChebyshev.okl", "SmoothChebyshevUpdate", kernelInfo);

    vectorAddInnerProdKernel = platform.buildKernel(PARALMOND_DIR"/okl/vectorAddInnerProd.okl", "vectorAddInnerProd", kernelInfo);

    kcycleCombinedOp1Kernel = platform.buildKernel(PARALMOND_DIR"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp1", kernelInfo);
    kcycleCombinedOp2Kernel = platform.buildKernel(PARALMOND_DIR"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp2", kernelInfo);

    dGEMVKernel = platform.buildKernel(PARALMOND_DIR"/okl/dGEMV.okl", "dGEMV", kernelInfo);

    if(rank==0) printf("done.\n");
  }
}

} //namespace parAlmond

} //namespace libp
