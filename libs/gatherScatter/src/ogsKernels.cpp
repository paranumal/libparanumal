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

#include "ogs.hpp"
#include "ogsKernels.hpp"

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
   //kernelInfo["compiler_flags"] += "-cl-opt-disable";
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
      ogs::gatherScatterKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatAdd", kernelInfo);
      ogs::gatherScatterKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMul", kernelInfo);
      ogs::gatherScatterKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMin", kernelInfo);
      ogs::gatherScatterKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_floatMax", kernelInfo);

      ogs::gatherScatterKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleAdd", kernelInfo);
      ogs::gatherScatterKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMul", kernelInfo);
      ogs::gatherScatterKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMin", kernelInfo);
      ogs::gatherScatterKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_doubleMax", kernelInfo);

      ogs::gatherScatterKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intAdd", kernelInfo);
      ogs::gatherScatterKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMul", kernelInfo);
      ogs::gatherScatterKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMin", kernelInfo);
      ogs::gatherScatterKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_intMax", kernelInfo);

      ogs::gatherScatterKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longAdd", kernelInfo);
      ogs::gatherScatterKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMul", kernelInfo);
      ogs::gatherScatterKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMin", kernelInfo);
      ogs::gatherScatterKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatter.okl", "gatherScatter_longMax", kernelInfo);

      ogs::gatherScatterVecKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatAdd", kernelInfo);
      ogs::gatherScatterVecKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMul", kernelInfo);
      ogs::gatherScatterVecKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMin", kernelInfo);
      ogs::gatherScatterVecKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_floatMax", kernelInfo);

      ogs::gatherScatterVecKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleAdd", kernelInfo);
      ogs::gatherScatterVecKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMul", kernelInfo);
      ogs::gatherScatterVecKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMin", kernelInfo);
      ogs::gatherScatterVecKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_doubleMax", kernelInfo);

      ogs::gatherScatterVecKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intAdd", kernelInfo);
      ogs::gatherScatterVecKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMul", kernelInfo);
      ogs::gatherScatterVecKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMin", kernelInfo);
      ogs::gatherScatterVecKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_intMax", kernelInfo);

      ogs::gatherScatterVecKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longAdd", kernelInfo);
      ogs::gatherScatterVecKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMul", kernelInfo);
      ogs::gatherScatterVecKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMin", kernelInfo);
      ogs::gatherScatterVecKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longMax", kernelInfo);

      ogs::gatherScatterManyKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatAdd", kernelInfo);
      ogs::gatherScatterManyKernel_floatMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMul", kernelInfo);
      ogs::gatherScatterManyKernel_floatMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMin", kernelInfo);
      ogs::gatherScatterManyKernel_floatMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_floatMax", kernelInfo);

      ogs::gatherScatterManyKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleAdd", kernelInfo);
      ogs::gatherScatterManyKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMul", kernelInfo);
      ogs::gatherScatterManyKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMin", kernelInfo);
      ogs::gatherScatterManyKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_doubleMax", kernelInfo);

      ogs::gatherScatterManyKernel_intAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intAdd", kernelInfo);
      ogs::gatherScatterManyKernel_intMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMul", kernelInfo);
      ogs::gatherScatterManyKernel_intMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMin", kernelInfo);
      ogs::gatherScatterManyKernel_intMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_intMax", kernelInfo);

      ogs::gatherScatterManyKernel_longAdd = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longAdd", kernelInfo);
      ogs::gatherScatterManyKernel_longMul = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMul", kernelInfo);
      ogs::gatherScatterManyKernel_longMin = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMin", kernelInfo);
      ogs::gatherScatterManyKernel_longMax = device.buildKernel(DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longMax", kernelInfo);
      


      ogs::gatherKernel_floatAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatAdd", kernelInfo);
      ogs::gatherKernel_floatMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMul", kernelInfo);
      ogs::gatherKernel_floatMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMin", kernelInfo);
      ogs::gatherKernel_floatMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_floatMax", kernelInfo);

      ogs::gatherKernel_doubleAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleAdd", kernelInfo);
      ogs::gatherKernel_doubleMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMul", kernelInfo);
      ogs::gatherKernel_doubleMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMin", kernelInfo);
      ogs::gatherKernel_doubleMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_doubleMax", kernelInfo);

      ogs::gatherKernel_intAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_intAdd", kernelInfo);
      ogs::gatherKernel_intMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMul", kernelInfo);
      ogs::gatherKernel_intMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMin", kernelInfo);
      ogs::gatherKernel_intMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_intMax", kernelInfo);

      ogs::gatherKernel_longAdd = device.buildKernel(DOGS "/okl/gather.okl", "gather_longAdd", kernelInfo);
      ogs::gatherKernel_longMul = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMul", kernelInfo);
      ogs::gatherKernel_longMin = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMin", kernelInfo);
      ogs::gatherKernel_longMax = device.buildKernel(DOGS "/okl/gather.okl", "gather_longMax", kernelInfo);

      ogs::gatherVecKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatAdd", kernelInfo);
      ogs::gatherVecKernel_floatMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMul", kernelInfo);
      ogs::gatherVecKernel_floatMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMin", kernelInfo);
      ogs::gatherVecKernel_floatMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_floatMax", kernelInfo);

      ogs::gatherVecKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleAdd", kernelInfo);
      ogs::gatherVecKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMul", kernelInfo);
      ogs::gatherVecKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMin", kernelInfo);
      ogs::gatherVecKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_doubleMax", kernelInfo);

      ogs::gatherVecKernel_intAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intAdd", kernelInfo);
      ogs::gatherVecKernel_intMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMul", kernelInfo);
      ogs::gatherVecKernel_intMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMin", kernelInfo);
      ogs::gatherVecKernel_intMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_intMax", kernelInfo);

      ogs::gatherVecKernel_longAdd = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longAdd", kernelInfo);
      ogs::gatherVecKernel_longMul = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMul", kernelInfo);
      ogs::gatherVecKernel_longMin = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMin", kernelInfo);
      ogs::gatherVecKernel_longMax = device.buildKernel(DOGS "/okl/gatherVec.okl", "gatherVec_longMax", kernelInfo);

      ogs::gatherManyKernel_floatAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatAdd", kernelInfo);
      ogs::gatherManyKernel_floatMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMul", kernelInfo);
      ogs::gatherManyKernel_floatMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMin", kernelInfo);
      ogs::gatherManyKernel_floatMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_floatMax", kernelInfo);

      ogs::gatherManyKernel_doubleAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleAdd", kernelInfo);
      ogs::gatherManyKernel_doubleMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMul", kernelInfo);
      ogs::gatherManyKernel_doubleMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMin", kernelInfo);
      ogs::gatherManyKernel_doubleMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_doubleMax", kernelInfo);

      ogs::gatherManyKernel_intAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intAdd", kernelInfo);
      ogs::gatherManyKernel_intMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMul", kernelInfo);
      ogs::gatherManyKernel_intMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMin", kernelInfo);
      ogs::gatherManyKernel_intMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_intMax", kernelInfo);

      ogs::gatherManyKernel_longAdd = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longAdd", kernelInfo);
      ogs::gatherManyKernel_longMul = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMul", kernelInfo);
      ogs::gatherManyKernel_longMin = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMin", kernelInfo);
      ogs::gatherManyKernel_longMax = device.buildKernel(DOGS "/okl/gatherMany.okl", "gatherMany_longMax", kernelInfo);

      

      ogs::scatterKernel_float = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_float", kernelInfo);
      ogs::scatterKernel_double = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_double", kernelInfo);
      ogs::scatterKernel_int = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_int", kernelInfo);
      ogs::scatterKernel_long = device.buildKernel(DOGS "/okl/scatter.okl", "scatter_long", kernelInfo);

      ogs::scatterVecKernel_float = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_float", kernelInfo);
      ogs::scatterVecKernel_double = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_double", kernelInfo);
      ogs::scatterVecKernel_int = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_int", kernelInfo);
      ogs::scatterVecKernel_long = device.buildKernel(DOGS "/okl/scatterVec.okl", "scatterVec_long", kernelInfo);

      ogs::scatterManyKernel_float = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_float", kernelInfo);
      ogs::scatterManyKernel_double = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_double", kernelInfo);
      ogs::scatterManyKernel_int = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_int", kernelInfo);
      ogs::scatterManyKernel_long = device.buildKernel(DOGS "/okl/scatterMany.okl", "scatterMany_long", kernelInfo);
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

