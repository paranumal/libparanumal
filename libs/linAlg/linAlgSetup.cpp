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

#include "core.hpp"
#include "linAlg.hpp"
#include "platform.hpp"

namespace libp {

linAlg_t::linAlg_t() {};

void linAlg_t::Setup(platform_t *_platform) {

  platform = _platform;
  kernelInfo = platform->props();

  //add defines
  kernelInfo["defines/" "p_blockSize"] = blocksize;
  kernelInfo["defines/init_dfloat_min"] =  std::numeric_limits<dfloat>::max();
  kernelInfo["defines/init_dfloat_max"] = -std::numeric_limits<dfloat>::max();

  //pinned scratch buffer
  h_scratch = platform->hostMalloc<dfloat>(blocksize);
  o_scratch = platform->malloc<dfloat>(blocksize);
}

//initialize list of kernels
void linAlg_t::InitKernels(std::vector<std::string> kernels) {

  for (size_t i=0;i<kernels.size();i++) {
    std::string name = kernels[i];
    if (name=="set") {
      if (setKernel.isInitialized()==false)
        setKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgSet.okl",
                                        "set",
                                        kernelInfo);
    } else if (name=="add") {
      if (addKernel.isInitialized()==false)
        addKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgAdd.okl",
                                        "add",
                                        kernelInfo);
    } else if (name=="scale") {
      if (scaleKernel.isInitialized()==false)
        scaleKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgScale.okl",
                                        "scale",
                                        kernelInfo);
    } else if (name=="axpy") {
      if (axpyKernel.isInitialized()==false)
        axpyKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgAXPY.okl",
                                        "axpy",
                                        kernelInfo);
    } else if (name=="zaxpy") {
      if (zaxpyKernel.isInitialized()==false)
        zaxpyKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgAXPY.okl",
                                        "zaxpy",
                                        kernelInfo);
    } else if (name=="amx") {
      if (amxKernel.isInitialized()==false)
        amxKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgAMXPY.okl",
                                        "amx",
                                        kernelInfo);
    } else if (name=="amxpy") {
      if (amxpyKernel.isInitialized()==false)
        amxpyKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgAMXPY.okl",
                                        "amxpy",
                                        kernelInfo);
    } else if (name=="zamxpy") {
      if (zamxpyKernel.isInitialized()==false)
        zamxpyKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgAMXPY.okl",
                                        "zamxpy",
                                        kernelInfo);
    } else if (name=="adx") {
      if (adxKernel.isInitialized()==false)
        adxKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgADXPY.okl",
                                        "adx",
                                        kernelInfo);
    } else if (name=="adxpy") {
      if (adxpyKernel.isInitialized()==false)
        adxpyKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgADXPY.okl",
                                        "adxpy",
                                        kernelInfo);
    } else if (name=="zadxpy") {
      if (zadxpyKernel.isInitialized()==false)
        zadxpyKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgADXPY.okl",
                                        "zadxpy",
                                        kernelInfo);
    } else if (name=="min") {
      if (minKernel1.isInitialized()==false) {
        minKernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgMin.okl",
                                        "min1",
                                        kernelInfo);
        minKernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgMin.okl",
                                        "min2",
                                        kernelInfo);
      }
    } else if (name=="max") {
      if (maxKernel1.isInitialized()==false) {
        maxKernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgMax.okl",
                                        "max1",
                                        kernelInfo);
        maxKernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgMax.okl",
                                        "max2",
                                        kernelInfo);
      }
    } else if (name=="sum") {
      if (sumKernel1.isInitialized()==false) {
        sumKernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgSum.okl",
                                        "sum1",
                                        kernelInfo);
        sumKernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgSum.okl",
                                        "sum2",
                                        kernelInfo);
      }
    } else if (name=="norm2") {
      if (norm2Kernel1.isInitialized()==false) {
        norm2Kernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgNorm2.okl",
                                        "norm2_1",
                                        kernelInfo);
        norm2Kernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgNorm2.okl",
                                        "norm2_2",
                                        kernelInfo);
      }
    } else if (name=="weightedNorm2") {
      if (weightedNorm2Kernel1.isInitialized()==false) {
        weightedNorm2Kernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgWeightedNorm2.okl",
                                        "weightedNorm2_1",
                                        kernelInfo);
        weightedNorm2Kernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgWeightedNorm2.okl",
                                        "weightedNorm2_1",
                                        kernelInfo);
      }
    } else if (name=="innerProd") {
      if (innerProdKernel1.isInitialized()==false) {
        innerProdKernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgInnerProd.okl",
                                        "innerProd1",
                                        kernelInfo);
        innerProdKernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgInnerProd.okl",
                                        "innerProd2",
                                        kernelInfo);
      }
    } else if (name=="weightedInnerProd") {
      if (weightedInnerProdKernel1.isInitialized()==false) {
        weightedInnerProdKernel1 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgWeightedInnerProd.okl",
                                        "weightedInnerProd1",
                                        kernelInfo);
        weightedInnerProdKernel2 = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgWeightedInnerProd.okl",
                                        "weightedInnerProd2",
                                        kernelInfo);
      }
    } else {
      LIBP_FORCE_ABORT("Requested linAlg routine \"" << name << "\" not found");
    }
  }
}

} //namespace libp
