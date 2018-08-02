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

#include "ogs.h"
#include "ogsKernels.h"

namespace ogs {

  int Nrefs = 0;

  void* hostBuf;
  size_t hostBufSize=0;

  void* haloBuf;
  occa::memory o_haloBuf;

  occa::kernel gatherScatterKernel;
  occa::kernel gatherKernel;
  occa::kernel scatterKernel;

  occa::stream defaultStream;
  occa::stream dataStream;

}


void ogs::initKernels(MPI_Comm comm, occa::device device) {

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  ogs::defaultStream = device.getStream();
  ogs::dataStream    = device.createStream();

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


  if(device.mode()=="OpenCL"){
   //    device.setCompilerFlags("-cl-opt-disable");
   kernelInfo["compiler_flags"] += "-cl-opt-disable";
  }

  if(device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
   kernelInfo["compiler_flags"] += "--ftz=true";
   kernelInfo["compiler_flags"] += "--prec-div=false";
   kernelInfo["compiler_flags"] += "--prec-sqrt=false";
   kernelInfo["compiler_flags"] += "--use_fast_math";
   kernelInfo["compiler_flags"] += "--fmad=true"; // compiler option for cuda
  }

  if (rank==0) printf("Compiling GatherScatter Kernels \n");

  for (int r=0;r<size;r++) {
    if (r==rank) {
      ogs::gatherScatterKernel 
        = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter", kernelInfo);
      
      ogs::gatherKernel 
        = device.buildKernel(DOGS "/okl/gather.okl", "gather", kernelInfo);

      ogs::scatterKernel 
        = device.buildKernel(DOGS "/okl/scatter.okl","scatter", kernelInfo);
    }
    MPI_Barrier(comm);
  }
}

void ogs::freeKernels() {

  ogs::gatherScatterKernel.free();
  ogs::gatherKernel.free();
  ogs::scatterKernel.free();

  ogs::o_haloBuf.free();
  ogs::haloBuf = NULL;
}

