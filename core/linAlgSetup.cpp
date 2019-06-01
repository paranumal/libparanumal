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

#include "linAlg.hpp"

#define LINALG_BLOCKSIZE 512

linAlg_t::linAlg_t(occa::device& device_,
         settings_t& settings_, occa::properties& props_):
  device(device_), settings(settings_), props(props_), blocksize(LINALG_BLOCKSIZE) {};

//named cosntructor
linAlg_t* linAlg_t::Setup(occa::device& device_,
         settings_t& settings_, occa::properties& props_) {

  linAlg_t *linAlg = new linAlg_t(device_, settings_, props_);

  //pinned scratch buffer
  occa::properties mprops;
  mprops["mapped"] = true;
  linAlg->h_scratch = linAlg->device.malloc(LINALG_BLOCKSIZE*sizeof(dfloat), mprops);
  linAlg->scratch = (dfloat*) linAlg->h_scratch.ptr();

  linAlg->o_scratch = linAlg->device.malloc(LINALG_BLOCKSIZE*sizeof(dfloat));

  return linAlg;
}

//initialize list of kernels
void linAlg_t::InitKernels(vector<string> kernels) {

  occa::properties kernelInfo = props; //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)LINALG_BLOCKSIZE;

  for (size_t i=0;i<kernels.size();i++) {
    string name = kernels[i];
    if (name=="set") {
      if (setKernel.isInitialized()==false)
        setKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgSet.okl",
                                        "set",
                                        kernelInfo);
    } else if (name=="scale") {
      if (scaleKernel.isInitialized()==false)
        scaleKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgScale.okl",
                                        "scale",
                                        kernelInfo);
    } else if (name=="add") {
      if (addKernel.isInitialized()==false)
        addKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAdd.okl",
                                        "add",
                                        kernelInfo);
    } else if (name=="axpy") {
      if (axpyKernel.isInitialized()==false)
        axpyKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAXPY.okl",
                                        "axpy",
                                        kernelInfo);
    } else if (name=="zaxpy") {
      if (zaxpyKernel.isInitialized()==false)
        zaxpyKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAXPY.okl",
                                        "zaxpy",
                                        kernelInfo);
    } else if (name=="axmy") {
      if (axmyKernel.isInitialized()==false)
        axmyKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAXMY.okl",
                                        "axmy",
                                        kernelInfo);
    } else if (name=="zaxmy") {
      if (zaxmyKernel.isInitialized()==false)
        zaxmyKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAXMY.okl",
                                        "zaxmy",
                                        kernelInfo);
    } else if (name=="axdy") {
      if (axdyKernel.isInitialized()==false)
        axdyKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAXDY.okl",
                                        "axdy",
                                        kernelInfo);
    } else if (name=="zaxdy") {
      if (zaxdyKernel.isInitialized()==false)
        zaxdyKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgAXDY.okl",
                                        "zaxdy",
                                        kernelInfo);
    } else if (name=="norm2") {
      if (norm2Kernel.isInitialized()==false)
        norm2Kernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgNorm2.okl",
                                        "norm2",
                                        kernelInfo);
    } else if (name=="innerProd") {
      if (innerProdKernel.isInitialized()==false)
        innerProdKernel = device.buildKernel(LIBP_DIR "/core/okl/"
                                        "linAlgInnerProd.okl",
                                        "innerProd",
                                        kernelInfo);
    } else if (name=="weightedInnerProd") {
      if (weightedInnerProdKernel.isInitialized()==false)
        weightedInnerProdKernel = device.buildKernel(LIBP_DIR "/core/okl/"
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
