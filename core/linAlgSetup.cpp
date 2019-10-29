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

#define LINALG_BLOCKSIZE 512

linAlg_t::linAlg_t(occa::device& device_,
         settings_t& settings_, occa::properties& props_):
  device(device_), settings(settings_), props(props_), blocksize(LINALG_BLOCKSIZE) {};

//named cosntructor
linAlg_t& linAlg_t::Setup(occa::device& device_,
         settings_t& settings_, occa::properties& props_) {

  linAlg_t *linAlg = new linAlg_t(device_, settings_, props_);

  //pinned scratch buffer
  occa::properties mprops;
  mprops["mapped"] = true;
  linAlg->h_scratch = linAlg->device.malloc(LINALG_BLOCKSIZE*sizeof(dfloat), mprops);
  linAlg->scratch = (dfloat*) linAlg->h_scratch.ptr(mprops);

  linAlg->o_scratch = linAlg->device.malloc(LINALG_BLOCKSIZE*sizeof(dfloat));

  return *linAlg;
}

//initialize list of kernels
void linAlg_t::InitKernels(vector<string> kernels, MPI_Comm& comm) {

  occa::properties kernelInfo = props; //copy base properties

  //add defines
  kernelInfo["defines/" "p_blockSize"] = (int)LINALG_BLOCKSIZE;

  for (size_t i=0;i<kernels.size();i++) {
    string name = kernels[i];
    if (name=="set") {
      if (setKernel.isInitialized()==false)
        setKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgSet.okl",
                                        "set",
                                        kernelInfo, comm);
    } else if (name=="add") {
      if (addKernel.isInitialized()==false)
        addKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgAdd.okl",
                                        "add",
                                        kernelInfo, comm);
    } else if (name=="scale") {
      if (scaleKernel.isInitialized()==false)
        scaleKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgScale.okl",
                                        "scale",
                                        kernelInfo, comm);
    } else if (name=="axpy") {
      if (axpyKernel.isInitialized()==false)
        axpyKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgAXPY.okl",
                                        "axpy",
                                        kernelInfo, comm);
    } else if (name=="zaxpy") {
      if (zaxpyKernel.isInitialized()==false)
        zaxpyKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgAXPY.okl",
                                        "zaxpy",
                                        kernelInfo, comm);
    } else if (name=="amx") {
      if (amxKernel.isInitialized()==false)
        amxKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgAMXPY.okl",
                                        "amx",
                                        kernelInfo, comm);
    } else if (name=="amxpy") {
      if (amxpyKernel.isInitialized()==false)
        amxpyKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgAMXPY.okl",
                                        "amxpy",
                                        kernelInfo, comm);
    } else if (name=="zamxpy") {
      if (zamxpyKernel.isInitialized()==false)
        zamxpyKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgAMXPY.okl",
                                        "zamxpy",
                                        kernelInfo, comm);
    } else if (name=="adx") {
      if (adxKernel.isInitialized()==false)
        adxKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgADXPY.okl",
                                        "adx",
                                        kernelInfo, comm);
    } else if (name=="adxpy") {
      if (adxpyKernel.isInitialized()==false)
        adxpyKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgADXPY.okl",
                                        "adxpy",
                                        kernelInfo, comm);
    } else if (name=="zadxpy") {
      if (zadxpyKernel.isInitialized()==false)
        zadxpyKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgADXPY.okl",
                                        "zadxpy",
                                        kernelInfo, comm);
    } else if (name=="sum") {
      if (sumKernel.isInitialized()==false)
        sumKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgSum.okl",
                                        "sum",
                                        kernelInfo, comm);
    } else if (name=="norm2") {
      if (norm2Kernel.isInitialized()==false)
        norm2Kernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgNorm2.okl",
                                        "norm2",
                                        kernelInfo, comm);
    } else if (name=="weightedNorm2") {
      if (weightedNorm2Kernel.isInitialized()==false)
        weightedNorm2Kernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgWeightedNorm2.okl",
                                        "weightedNorm2",
                                        kernelInfo, comm);
    } else if (name=="innerProd") {
      if (innerProdKernel.isInitialized()==false)
        innerProdKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgInnerProd.okl",
                                        "innerProd",
                                        kernelInfo, comm);
    } else if (name=="weightedInnerProd") {
      if (weightedInnerProdKernel.isInitialized()==false)
        weightedInnerProdKernel = buildKernel(device, LIBP_DIR "/core/okl/"
                                        "linAlgWeightedInnerProd.okl",
                                        "weightedInnerProd",
                                        kernelInfo, comm);
    } else {
      stringstream ss;
      ss << "Requested linAlg routine \"" << name << "\" not found";
      LIBP_ABORT(ss.str());
    }
  }
}

linAlg_t::~linAlg_t() {
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
  sumKernel.free();
  norm2Kernel.free();
  weightedNorm2Kernel.free();
  innerProdKernel.free();
  weightedInnerProdKernel.free();
}