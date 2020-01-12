/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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
#include "parAlmond.hpp"

namespace parAlmond {

int Nrefs = 0;

occa::kernel haloExtractKernel;

occa::kernel SpMVcsrKernel1;
occa::kernel SpMVcsrKernel2;
occa::kernel SpMVellKernel1;
occa::kernel SpMVellKernel2;
occa::kernel SpMVmcsrKernel1;
occa::kernel SpMVmcsrKernel2;

occa::kernel vectorSetKernel;
occa::kernel vectorScaleKernel;
occa::kernel vectorAddScalarKernel;
occa::kernel vectorAddKernel1;
occa::kernel vectorAddKernel2;
occa::kernel vectorDotStarKernel1;
occa::kernel vectorDotStarKernel2;
occa::kernel vectorInnerProdKernel;
occa::kernel kcycleCombinedOp1Kernel;
occa::kernel kcycleCombinedOp2Kernel;
occa::kernel kcycleWeightedCombinedOp1Kernel;
occa::kernel kcycleWeightedCombinedOp2Kernel;
occa::kernel vectorAddInnerProdKernel;
occa::kernel vectorAddWeightedInnerProdKernel;

void buildParAlmondKernels(MPI_Comm comm, occa::device device){

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  double seed = (double) rank;
  srand48(seed);

  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  if(sizeof(dlong)==4){
    kernelInfo["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
    kernelInfo["defines/" "dlong"]="long long int";
  }

  if(sizeof(dfloat) == sizeof(double)){
    kernelInfo["defines/" "dfloat"]= "double";
    kernelInfo["defines/" "dfloat4"]= "double4";
  }
  else if(sizeof(dfloat) == sizeof(float)){
    kernelInfo["defines/" "dfloat"]= "float";
    kernelInfo["defines/" "dfloat4"]= "float4";
  }

  kernelInfo["defines/" "p_BLOCKSIZE"]= BLOCKSIZE;

  if(device.mode()=="OpenCL"){
    //kernelInfo["compiler_flags"] += "-cl-opt-disable";
  }

  if(device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    kernelInfo["compiler_flags"] += " --ftz=true";
    kernelInfo["compiler_flags"] += " --prec-div=false";
    kernelInfo["compiler_flags"] += " --prec-sqrt=false";
    kernelInfo["compiler_flags"] += " --use_fast_math";
    kernelInfo["compiler_flags"] += " --fmad=true"; // compiler option for cuda
  }

  if (rank==0) {printf("Compiling parALMOND Kernels...");fflush(stdout);}

  SpMVcsrKernel1  = buildKernel(device, DPARALMOND"/okl/SpMVcsr.okl",  "SpMVcsr1",  kernelInfo, comm);
  SpMVcsrKernel2  = buildKernel(device, DPARALMOND"/okl/SpMVcsr.okl",  "SpMVcsr2",  kernelInfo, comm);
  SpMVellKernel1  = buildKernel(device, DPARALMOND"/okl/SpMVell.okl",  "SpMVell1",  kernelInfo, comm);
  SpMVellKernel2  = buildKernel(device, DPARALMOND"/okl/SpMVell.okl",  "SpMVell2",  kernelInfo, comm);
  SpMVmcsrKernel1 = buildKernel(device, DPARALMOND"/okl/SpMVmcsr.okl", "SpMVmcsr1", kernelInfo, comm);
  SpMVmcsrKernel2 = buildKernel(device, DPARALMOND"/okl/SpMVmcsr.okl", "SpMVmcsr2", kernelInfo, comm);

  vectorSetKernel = buildKernel(device, DPARALMOND"/okl/vectorSet.okl", "vectorSet", kernelInfo, comm);
  vectorScaleKernel = buildKernel(device, DPARALMOND"/okl/vectorScale.okl", "vectorScale", kernelInfo, comm);
  vectorAddScalarKernel = buildKernel(device, DPARALMOND"/okl/vectorAddScalar.okl", "vectorAddScalar", kernelInfo, comm);
  vectorAddKernel1 = buildKernel(device, DPARALMOND"/okl/vectorAdd.okl", "vectorAdd1", kernelInfo, comm);
  vectorAddKernel2 = buildKernel(device, DPARALMOND"/okl/vectorAdd.okl", "vectorAdd2", kernelInfo, comm);
  vectorDotStarKernel1 = buildKernel(device, DPARALMOND"/okl/vectorDotStar.okl", "vectorDotStar1", kernelInfo, comm);
  vectorDotStarKernel2 = buildKernel(device, DPARALMOND"/okl/vectorDotStar.okl", "vectorDotStar2", kernelInfo, comm);
  vectorInnerProdKernel = buildKernel(device, DPARALMOND"/okl/vectorInnerProd.okl", "vectorInnerProd", kernelInfo, comm);

  vectorAddInnerProdKernel = buildKernel(device, DPARALMOND"/okl/vectorAddInnerProd.okl", "vectorAddInnerProd", kernelInfo, comm);
  vectorAddWeightedInnerProdKernel = buildKernel(device, DPARALMOND"/okl/vectorAddInnerProd.okl", "vectorAddWeightedInnerProd", kernelInfo, comm);

  kcycleCombinedOp1Kernel = buildKernel(device, DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp1", kernelInfo, comm);
  kcycleCombinedOp2Kernel = buildKernel(device, DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleCombinedOp2", kernelInfo, comm);
  kcycleWeightedCombinedOp1Kernel = buildKernel(device, DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleWeightedCombinedOp1", kernelInfo, comm);
  kcycleWeightedCombinedOp2Kernel = buildKernel(device, DPARALMOND"/okl/kcycleCombinedOp.okl", "kcycleWeightedCombinedOp2", kernelInfo, comm);

  haloExtractKernel = buildKernel(device, DPARALMOND"/okl/haloExtract.okl", "haloExtract", kernelInfo, comm);

  if(rank==0) printf("done.\n");
}

void freeParAlmondKernels() {

  haloExtractKernel.free();

  SpMVcsrKernel1.free();
  SpMVcsrKernel2.free();
  SpMVellKernel1.free();
  SpMVellKernel2.free();
  SpMVmcsrKernel1.free();
  SpMVmcsrKernel2.free();

  vectorSetKernel.free();
  vectorScaleKernel.free();
  vectorAddScalarKernel.free();
  vectorAddKernel1.free();
  vectorAddKernel2.free();
  vectorDotStarKernel1.free();
  vectorDotStarKernel2.free();
  vectorInnerProdKernel.free();
  kcycleCombinedOp1Kernel.free();
  kcycleCombinedOp2Kernel.free();
  kcycleWeightedCombinedOp1Kernel.free();
  kcycleWeightedCombinedOp2Kernel.free();
  vectorAddInnerProdKernel.free();
  vectorAddWeightedInnerProdKernel.free();

}


} //namespace parAlmond
