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

#include "core.hpp"
#include "linAlg.hpp"
#include "platform.hpp"

#define LINALG_BLOCKSIZE 512

template <class lafloat>
linAlg_t<lafloat>::linAlg_t(): blocksize(LINALG_BLOCKSIZE) {};

template <class lafloat>
void linAlg_t<lafloat>::Setup(platform_t *_platform) {

  platform = _platform;
  kernelInfo = platform->props;

  if(sizeof(lafloat)==sizeof(float)){
    printf("build linAlg kernels with floats\n");
    kernelInfo["defines/dfloat"] = "float";
  }
  if(sizeof(lafloat)==sizeof(double)){
    kernelInfo["defines/dfloat"] = "double";
  }
  
  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)LINALG_BLOCKSIZE;

  kernelInfo["defines/init_dfloat_min"] =  std::numeric_limits<dfloat>::max();
  kernelInfo["defines/init_dfloat_max"] = -std::numeric_limits<dfloat>::max();

  //pinned scratch buffer
  scratch = (lafloat*) platform->hostMalloc(LINALG_BLOCKSIZE*sizeof(lafloat),
                                           NULL, h_scratch);
  o_scratch = platform->malloc(LINALG_BLOCKSIZE*sizeof(lafloat));
}

//initialize list of kernels
template <class lafloat>
void linAlg_t<lafloat>::InitKernels(vector<string> kernels) {

  for (size_t i=0;i<kernels.size();i++) {
    string name = kernels[i];
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
      if (minKernel.isInitialized()==false)
        minKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgMin.okl",
                                        "min",
                                        kernelInfo);
    } else if (name=="max") {
      if (maxKernel.isInitialized()==false)
        maxKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgMax.okl",
                                        "max",
                                        kernelInfo);
    } else if (name=="sum") {
      if (sumKernel.isInitialized()==false)
        sumKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgSum.okl",
                                        "sum",
                                        kernelInfo);
    } else if (name=="norm2") {
      if (norm2Kernel.isInitialized()==false)
        norm2Kernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgNorm2.okl",
                                        "norm2",
                                        kernelInfo);
    } else if (name=="weightedNorm2") {
      if (weightedNorm2Kernel.isInitialized()==false)
        weightedNorm2Kernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgWeightedNorm2.okl",
                                        "weightedNorm2",
                                        kernelInfo);
    } else if (name=="innerProd") {
      if (innerProdKernel.isInitialized()==false)
        innerProdKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgInnerProd.okl",
                                        "innerProd",
                                        kernelInfo);
    } else if (name=="weightedInnerProd") {
      if (weightedInnerProdKernel.isInitialized()==false)
        weightedInnerProdKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                        "linAlgWeightedInnerProd.okl",
                                        "weightedInnerProd",
                                        kernelInfo);
    } else {
      stringstream ss;
      ss << "Requested linAlg routine \"" << name << "\" not found";
      LIBP_ABORT(ss.str());
    }
  }
}

template <class lafloat>
linAlg_t<lafloat>::~linAlg_t() {
  setKernel.free();
  addKernel.free();
  scaleKernel.free();
  axpyKernel.free();
  zaxpyKernel.free();
  amxKernel.free();
  amxpyKernel.free();
  zamxpyKernel.free();
  adxKernel.free();
  adxpyKernel.free();
  zadxpyKernel.free();
  minKernel.free();
  maxKernel.free();
  sumKernel.free();
  norm2Kernel.free();
  weightedNorm2Kernel.free();
  innerProdKernel.free();
  weightedInnerProdKernel.free();
}

template class linAlg_t<float>;
template class linAlg_t<double>;

//template class void linAlg_t<float>::Setup(platform_t *_platform);
//template class void linAlg_t<double>::Setup(platform_t *_platform);
