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
  kernelInfoFloat = platform->props();

  //add defines
  kernelInfoFloat["defines/" "p_blockSize"] = blocksize;
  kernelInfoFloat["defines/init_dfloat_min"] =  std::numeric_limits<float>::max();
  kernelInfoFloat["defines/init_dfloat_max"] = -std::numeric_limits<float>::max();
  kernelInfoFloat["defines/dfloat"]  = "float";

  kernelInfoDouble = kernelInfoFloat;
  kernelInfoDouble["defines/init_dfloat_min"] =  std::numeric_limits<double>::max();
  kernelInfoDouble["defines/init_dfloat_max"] = -std::numeric_limits<double>::max();
  kernelInfoDouble["defines/dfloat"]  = "double";
}

//initialize list of kernels
void linAlg_t::InitKernels(std::vector<std::string> kernels) {

  for (size_t i=0;i<kernels.size();i++) {
    std::string name = kernels[i];
    if (name=="set") {
      if (setKernelFloat.isInitialized()==false){
        setKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgSet.okl",
                                               "set",
                                               kernelInfoFloat);
        setKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgSet.okl",
                                                "set",
                                                kernelInfoDouble);

      }
    } else if (name=="add") {
      if (addKernelFloat.isInitialized()==false){
        addKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgAdd.okl",
                                               "add",
                                               kernelInfoFloat);
        addKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgAdd.okl",
                                                "add",
                                                kernelInfoDouble);
      }
    } else if (name=="scale") {
      if (scaleKernelFloat.isInitialized()==false){
        scaleKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                 "linAlgScale.okl",
                                                 "scale",
                                                 kernelInfoFloat);
        scaleKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgScale.okl",
                                                  "scale",
                                                  kernelInfoDouble);
      }

    } else if (name=="axpy") {
      if (axpyKernelFloat.isInitialized()==false){
        axpyKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgAXPY.okl",
                                                "axpy",
                                                kernelInfoFloat);
        axpyKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                 "linAlgAXPY.okl",
                                                 "axpy",
                                                 kernelInfoDouble);
      }
    } else if (name=="zaxpy") {
      if (zaxpyKernelFloat.isInitialized()==false){
        zaxpyKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                 "linAlgAXPY.okl",
                                                 "zaxpy",
                                                 kernelInfoFloat);
        zaxpyKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgAXPY.okl",
                                                  "zaxpy",
                                                  kernelInfoDouble);
      }
    } else if (name=="amx") {
      if (amxKernelFloat.isInitialized()==false){
        amxKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgAMXPY.okl",
                                               "amx",
                                               kernelInfoFloat);

        amxKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgAMXPY.okl",
                                                "amx",
                                                kernelInfoDouble);
      }

    } else if (name=="amxpy") {
      if (amxpyKernelFloat.isInitialized()==false){
        amxpyKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                 "linAlgAMXPY.okl",
                                                 "amxpy",
                                                 kernelInfoFloat);
        amxpyKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgAMXPY.okl",
                                                  "amxpy",
                                                  kernelInfoDouble);
      }
    } else if (name=="zamxpy") {
      if (zamxpyKernelFloat.isInitialized()==false){
        zamxpyKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgAMXPY.okl",
                                                  "zamxpy",
                                                  kernelInfoFloat);
        zamxpyKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                   "linAlgAMXPY.okl",
                                                   "zamxpy",
                                                   kernelInfoDouble);
      }
    } else if (name=="adx") {
      if (adxKernelFloat.isInitialized()==false){
        adxKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgADXPY.okl",
                                               "adx",
                                               kernelInfoFloat);
        adxKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgADXPY.okl",
                                                "adx",
                                                kernelInfoDouble);
      }
    } else if (name=="adxpy") {
      if (adxpyKernelFloat.isInitialized()==false){
        adxpyKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                 "linAlgADXPY.okl",
                                                 "adxpy",
                                                 kernelInfoFloat);
        adxpyKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgADXPY.okl",
                                                  "adxpy",
                                                  kernelInfoDouble);
      }
    } else if (name=="zadxpy") {
      if (zadxpyKernelFloat.isInitialized()==false){
        zadxpyKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgADXPY.okl",
                                                  "zadxpy",
                                                  kernelInfoFloat);
        zadxpyKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                   "linAlgADXPY.okl",
                                                   "zadxpy",
                                                   kernelInfoDouble);
      }
    } else if (name=="min") {
      if (minKernelFloat.isInitialized()==false) {
        minKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgMin.okl",
                                               "min",
                                               kernelInfoFloat);
        minKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgMin.okl",
                                                "min",
                                                kernelInfoDouble);
      }
    } else if (name=="max") {
      if (maxKernelFloat.isInitialized()==false) {
        maxKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgMax.okl",
                                               "max",
                                               kernelInfoFloat);
        maxKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgMax.okl",
                                                "max",
                                                kernelInfoDouble);
      }
    } else if (name=="sum") {
      if (sumKernelFloat.isInitialized()==false) {
        sumKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                               "linAlgSum.okl",
                                               "sum",
                                               kernelInfoFloat);
        sumKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                "linAlgSum.okl",
                                                "sum",
                                                kernelInfoDouble);
      }
    } else if (name=="norm2") {
      if (norm2KernelFloat.isInitialized()==false) {
        norm2KernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                 "linAlgNorm2.okl",
                                                 "norm2",
                                                 kernelInfoFloat);
        norm2KernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                  "linAlgNorm2.okl",
                                                  "norm2",
                                                  kernelInfoDouble);
      }
    } else if (name=="weightedNorm2") {
      if (weightedNorm2KernelFloat.isInitialized()==false) {
        weightedNorm2KernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                         "linAlgWeightedNorm2.okl",
                                                         "weightedNorm2",
                                                         kernelInfoFloat);
        weightedNorm2KernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                          "linAlgWeightedNorm2.okl",
                                                          "weightedNorm2",
                                                          kernelInfoDouble);
      }
    } else if (name=="innerProd") {
      if (innerProdKernelFloat.isInitialized()==false) {
        innerProdKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                     "linAlgInnerProd.okl",
                                                     "innerProd",
                                                     kernelInfoFloat);
        innerProdKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                      "linAlgInnerProd.okl",
                                                      "innerProd",
                                                      kernelInfoDouble);
      }
    } else if (name=="weightedInnerProd") {
      if (weightedInnerProdKernelFloat.isInitialized()==false) {
        weightedInnerProdKernelFloat = platform->buildKernel(LINALG_DIR "/okl/"
                                                             "linAlgWeightedInnerProd.okl",
                                                             "weightedInnerProd",
                                                             kernelInfoFloat);
        weightedInnerProdKernelDouble = platform->buildKernel(LINALG_DIR "/okl/"
                                                              "linAlgWeightedInnerProd.okl",
                                                              "weightedInnerProd",
                                                              kernelInfoDouble);
      }
    } else if (name=="d2p" || name=="p2d"){
      properties_t kernelInfoType = platform->props();

      //add defines
      kernelInfoType["defines/" "p_blockSize"] = blocksize;
      kernelInfoType["defines/dfloat"]  = dfloatString;
      kernelInfoType["defines/pfloat"]  = pfloatString;

      if(d2pKernel.isInitialized()==false){
        d2pKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                          "linAlgTypeConvert.okl",
                                          "d2p",
                                          kernelInfoType);
      }

      if(p2dKernel.isInitialized()==false){
        p2dKernel = platform->buildKernel(LINALG_DIR "/okl/"
                                          "linAlgTypeConvert.okl",
                                          "p2d",
                                          kernelInfoType);
      }
    } else {
      LIBP_FORCE_ABORT("Requested linAlg routine \"" << name << "\" not found");
    }
  }
}

} //namespace libp
