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
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include <limits>

namespace ogs {

  int Nrefs = 0;

  void* hostBuf;
  size_t hostBufSize=0;

  void* haloBuf;
  occa::memory o_haloBuf;
  occa::memory h_haloBuf;

  occa::stream defaultStream;
  occa::stream dataStream;

  occa::kernel gatherScatterKernel_floatAdd;
  occa::kernel gatherScatterKernel_floatMul;
  occa::kernel gatherScatterKernel_floatMin;
  occa::kernel gatherScatterKernel_floatMax;
  occa::kernel gatherScatterKernel_doubleAdd;
  occa::kernel gatherScatterKernel_doubleMul;
  occa::kernel gatherScatterKernel_doubleMin;
  occa::kernel gatherScatterKernel_doubleMax;
  occa::kernel gatherScatterKernel_intAdd;
  occa::kernel gatherScatterKernel_intMul;
  occa::kernel gatherScatterKernel_intMin;
  occa::kernel gatherScatterKernel_intMax;
  occa::kernel gatherScatterKernel_longAdd;
  occa::kernel gatherScatterKernel_longMul;
  occa::kernel gatherScatterKernel_longMin;
  occa::kernel gatherScatterKernel_longMax;
  occa::kernel gatherScatterVecKernel_floatAdd;
  occa::kernel gatherScatterVecKernel_floatMul;
  occa::kernel gatherScatterVecKernel_floatMin;
  occa::kernel gatherScatterVecKernel_floatMax;
  occa::kernel gatherScatterVecKernel_doubleAdd;
  occa::kernel gatherScatterVecKernel_doubleMul;
  occa::kernel gatherScatterVecKernel_doubleMin;
  occa::kernel gatherScatterVecKernel_doubleMax;
  occa::kernel gatherScatterVecKernel_intAdd;
  occa::kernel gatherScatterVecKernel_intMul;
  occa::kernel gatherScatterVecKernel_intMin;
  occa::kernel gatherScatterVecKernel_intMax;
  occa::kernel gatherScatterVecKernel_longAdd;
  occa::kernel gatherScatterVecKernel_longMul;
  occa::kernel gatherScatterVecKernel_longMin;
  occa::kernel gatherScatterVecKernel_longMax;
  occa::kernel gatherScatterManyKernel_floatAdd;
  occa::kernel gatherScatterManyKernel_floatMul;
  occa::kernel gatherScatterManyKernel_floatMin;
  occa::kernel gatherScatterManyKernel_floatMax;
  occa::kernel gatherScatterManyKernel_doubleAdd;
  occa::kernel gatherScatterManyKernel_doubleMul;
  occa::kernel gatherScatterManyKernel_doubleMin;
  occa::kernel gatherScatterManyKernel_doubleMax;
  occa::kernel gatherScatterManyKernel_intAdd;
  occa::kernel gatherScatterManyKernel_intMul;
  occa::kernel gatherScatterManyKernel_intMin;
  occa::kernel gatherScatterManyKernel_intMax;
  occa::kernel gatherScatterManyKernel_longAdd;
  occa::kernel gatherScatterManyKernel_longMul;
  occa::kernel gatherScatterManyKernel_longMin;
  occa::kernel gatherScatterManyKernel_longMax;

  occa::kernel gatherKernel_floatAdd;
  occa::kernel gatherKernel_floatMul;
  occa::kernel gatherKernel_floatMin;
  occa::kernel gatherKernel_floatMax;
  occa::kernel gatherKernel_doubleAdd;
  occa::kernel gatherKernel_doubleMul;
  occa::kernel gatherKernel_doubleMin;
  occa::kernel gatherKernel_doubleMax;
  occa::kernel gatherKernel_intAdd;
  occa::kernel gatherKernel_intMul;
  occa::kernel gatherKernel_intMin;
  occa::kernel gatherKernel_intMax;
  occa::kernel gatherKernel_longAdd;
  occa::kernel gatherKernel_longMul;
  occa::kernel gatherKernel_longMin;
  occa::kernel gatherKernel_longMax;
  occa::kernel gatherVecKernel_floatAdd;
  occa::kernel gatherVecKernel_floatMul;
  occa::kernel gatherVecKernel_floatMin;
  occa::kernel gatherVecKernel_floatMax;
  occa::kernel gatherVecKernel_doubleAdd;
  occa::kernel gatherVecKernel_doubleMul;
  occa::kernel gatherVecKernel_doubleMin;
  occa::kernel gatherVecKernel_doubleMax;
  occa::kernel gatherVecKernel_intAdd;
  occa::kernel gatherVecKernel_intMul;
  occa::kernel gatherVecKernel_intMin;
  occa::kernel gatherVecKernel_intMax;
  occa::kernel gatherVecKernel_longAdd;
  occa::kernel gatherVecKernel_longMul;
  occa::kernel gatherVecKernel_longMin;
  occa::kernel gatherVecKernel_longMax;
  occa::kernel gatherManyKernel_floatAdd;
  occa::kernel gatherManyKernel_floatMul;
  occa::kernel gatherManyKernel_floatMin;
  occa::kernel gatherManyKernel_floatMax;
  occa::kernel gatherManyKernel_doubleAdd;
  occa::kernel gatherManyKernel_doubleMul;
  occa::kernel gatherManyKernel_doubleMin;
  occa::kernel gatherManyKernel_doubleMax;
  occa::kernel gatherManyKernel_intAdd;
  occa::kernel gatherManyKernel_intMul;
  occa::kernel gatherManyKernel_intMin;
  occa::kernel gatherManyKernel_intMax;
  occa::kernel gatherManyKernel_longAdd;
  occa::kernel gatherManyKernel_longMul;
  occa::kernel gatherManyKernel_longMin;
  occa::kernel gatherManyKernel_longMax;

  occa::kernel scatterKernel_float;
  occa::kernel scatterKernel_double;
  occa::kernel scatterKernel_int;
  occa::kernel scatterKernel_long;
  occa::kernel scatterVecKernel_float;
  occa::kernel scatterVecKernel_double;
  occa::kernel scatterVecKernel_int;
  occa::kernel scatterVecKernel_long;
  occa::kernel scatterManyKernel_float;
  occa::kernel scatterManyKernel_double;
  occa::kernel scatterManyKernel_int;
  occa::kernel scatterManyKernel_long;
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

  kernelInfo["defines/init_float_add"] = (float) 0.0;
  kernelInfo["defines/init_float_mul"] = (float) 1.0;
  kernelInfo["defines/init_float_min"] = (float)  std::numeric_limits<float>::max();
  kernelInfo["defines/init_float_max"] = (float) -std::numeric_limits<float>::max();

  kernelInfo["defines/init_double_add"] = (double) 0.0;
  kernelInfo["defines/init_double_mul"] = (double) 1.0;
  kernelInfo["defines/init_double_min"] = (double)  std::numeric_limits<double>::max();
  kernelInfo["defines/init_double_max"] = (double) -std::numeric_limits<double>::max();

  kernelInfo["defines/init_int_add"] = (int) 0;
  kernelInfo["defines/init_int_mul"] = (int) 1;
  kernelInfo["defines/init_int_min"] = (int)  std::numeric_limits<int>::max();
  kernelInfo["defines/init_int_max"] = (int) -std::numeric_limits<int>::max();

  kernelInfo["defines/init_long_add"] = (long) 0;
  kernelInfo["defines/init_long_mul"] = (long) 1;
  kernelInfo["defines/init_long_min"] = (long)  std::numeric_limits<long>::max();
  kernelInfo["defines/init_long_max"] = (long) -std::numeric_limits<long>::max();

  kernelInfo["defines/init_longlongint_add"] = (int64_t) 0;
  kernelInfo["defines/init_longlongint_mul"] = (int64_t) 1;
  kernelInfo["defines/init_longlongint_min"] = (int64_t)  std::numeric_limits<int64_t>::max();
  kernelInfo["defines/init_longlongint_max"] = (int64_t) -std::numeric_limits<int64_t>::max();

  kernelInfo["includes"] += DOGS "/okl/ogsDefs.h";

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

  if (rank==0) {printf("Compiling GatherScatter Kernels...");fflush(stdout);}

  ogs::gatherScatterKernel_floatAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_float_add", kernelInfo, comm);
  ogs::gatherScatterKernel_floatMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_float_mul", kernelInfo, comm);
  ogs::gatherScatterKernel_floatMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_float_min", kernelInfo, comm);
  ogs::gatherScatterKernel_floatMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_float_max", kernelInfo, comm);

  ogs::gatherScatterKernel_doubleAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_double_add", kernelInfo, comm);
  ogs::gatherScatterKernel_doubleMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_double_mul", kernelInfo, comm);
  ogs::gatherScatterKernel_doubleMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_double_min", kernelInfo, comm);
  ogs::gatherScatterKernel_doubleMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_double_max", kernelInfo, comm);

  ogs::gatherScatterKernel_intAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_int_add", kernelInfo, comm);
  ogs::gatherScatterKernel_intMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_int_mul", kernelInfo, comm);
  ogs::gatherScatterKernel_intMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_int_min", kernelInfo, comm);
  ogs::gatherScatterKernel_intMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_int_max", kernelInfo, comm);

  ogs::gatherScatterKernel_longAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_longlongint_add", kernelInfo, comm);
  ogs::gatherScatterKernel_longMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_longlongint_mul", kernelInfo, comm);
  ogs::gatherScatterKernel_longMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_longlongint_min", kernelInfo, comm);
  ogs::gatherScatterKernel_longMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gatherScatter_longlongint_max", kernelInfo, comm);

  ogs::gatherScatterVecKernel_floatAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_float_add", kernelInfo, comm);
  ogs::gatherScatterVecKernel_floatMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_float_mul", kernelInfo, comm);
  ogs::gatherScatterVecKernel_floatMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_float_min", kernelInfo, comm);
  ogs::gatherScatterVecKernel_floatMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_float_max", kernelInfo, comm);

  ogs::gatherScatterVecKernel_doubleAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_double_add", kernelInfo, comm);
  ogs::gatherScatterVecKernel_doubleMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_double_mul", kernelInfo, comm);
  ogs::gatherScatterVecKernel_doubleMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_double_min", kernelInfo, comm);
  ogs::gatherScatterVecKernel_doubleMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_double_max", kernelInfo, comm);

  ogs::gatherScatterVecKernel_intAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_int_add", kernelInfo, comm);
  ogs::gatherScatterVecKernel_intMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_int_mul", kernelInfo, comm);
  ogs::gatherScatterVecKernel_intMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_int_min", kernelInfo, comm);
  ogs::gatherScatterVecKernel_intMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_int_max", kernelInfo, comm);

  ogs::gatherScatterVecKernel_longAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longlongint_add", kernelInfo, comm);
  ogs::gatherScatterVecKernel_longMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longlongint_mul", kernelInfo, comm);
  ogs::gatherScatterVecKernel_longMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longlongint_min", kernelInfo, comm);
  ogs::gatherScatterVecKernel_longMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherScatterVec_longlongint_max", kernelInfo, comm);

  ogs::gatherScatterManyKernel_floatAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_float_add", kernelInfo, comm);
  ogs::gatherScatterManyKernel_floatMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_float_mul", kernelInfo, comm);
  ogs::gatherScatterManyKernel_floatMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_float_min", kernelInfo, comm);
  ogs::gatherScatterManyKernel_floatMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_float_max", kernelInfo, comm);

  ogs::gatherScatterManyKernel_doubleAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_double_add", kernelInfo, comm);
  ogs::gatherScatterManyKernel_doubleMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_double_mul", kernelInfo, comm);
  ogs::gatherScatterManyKernel_doubleMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_double_min", kernelInfo, comm);
  ogs::gatherScatterManyKernel_doubleMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_double_max", kernelInfo, comm);

  ogs::gatherScatterManyKernel_intAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_int_add", kernelInfo, comm);
  ogs::gatherScatterManyKernel_intMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_int_mul", kernelInfo, comm);
  ogs::gatherScatterManyKernel_intMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_int_min", kernelInfo, comm);
  ogs::gatherScatterManyKernel_intMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_int_max", kernelInfo, comm);

  ogs::gatherScatterManyKernel_longAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longlongint_add", kernelInfo, comm);
  ogs::gatherScatterManyKernel_longMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longlongint_mul", kernelInfo, comm);
  ogs::gatherScatterManyKernel_longMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longlongint_min", kernelInfo, comm);
  ogs::gatherScatterManyKernel_longMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherScatterMany_longlongint_max", kernelInfo, comm);



  ogs::gatherKernel_floatAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_float_add", kernelInfo, comm);
  ogs::gatherKernel_floatMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_float_mul", kernelInfo, comm);
  ogs::gatherKernel_floatMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_float_min", kernelInfo, comm);
  ogs::gatherKernel_floatMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_float_max", kernelInfo, comm);

  ogs::gatherKernel_doubleAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_double_add", kernelInfo, comm);
  ogs::gatherKernel_doubleMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_double_mul", kernelInfo, comm);
  ogs::gatherKernel_doubleMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_double_min", kernelInfo, comm);
  ogs::gatherKernel_doubleMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_double_max", kernelInfo, comm);

  ogs::gatherKernel_intAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_int_add", kernelInfo, comm);
  ogs::gatherKernel_intMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_int_mul", kernelInfo, comm);
  ogs::gatherKernel_intMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_int_min", kernelInfo, comm);
  ogs::gatherKernel_intMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_int_max", kernelInfo, comm);

  ogs::gatherKernel_longAdd = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_longlongint_add", kernelInfo, comm);
  ogs::gatherKernel_longMul = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_longlongint_mul", kernelInfo, comm);
  ogs::gatherKernel_longMin = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_longlongint_min", kernelInfo, comm);
  ogs::gatherKernel_longMax = buildKernel(device, DOGS "/okl/gatherScatter.okl", "gather_longlongint_max", kernelInfo, comm);

  ogs::gatherVecKernel_floatAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_float_add", kernelInfo, comm);
  ogs::gatherVecKernel_floatMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_float_mul", kernelInfo, comm);
  ogs::gatherVecKernel_floatMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_float_min", kernelInfo, comm);
  ogs::gatherVecKernel_floatMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_float_max", kernelInfo, comm);

  ogs::gatherVecKernel_doubleAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_double_add", kernelInfo, comm);
  ogs::gatherVecKernel_doubleMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_double_mul", kernelInfo, comm);
  ogs::gatherVecKernel_doubleMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_double_min", kernelInfo, comm);
  ogs::gatherVecKernel_doubleMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_double_max", kernelInfo, comm);

  ogs::gatherVecKernel_intAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_int_add", kernelInfo, comm);
  ogs::gatherVecKernel_intMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_int_mul", kernelInfo, comm);
  ogs::gatherVecKernel_intMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_int_min", kernelInfo, comm);
  ogs::gatherVecKernel_intMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_int_max", kernelInfo, comm);

  ogs::gatherVecKernel_longAdd = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_longlongint_add", kernelInfo, comm);
  ogs::gatherVecKernel_longMul = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_longlongint_mul", kernelInfo, comm);
  ogs::gatherVecKernel_longMin = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_longlongint_min", kernelInfo, comm);
  ogs::gatherVecKernel_longMax = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "gatherVec_longlongint_max", kernelInfo, comm);

  ogs::gatherManyKernel_floatAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_float_add", kernelInfo, comm);
  ogs::gatherManyKernel_floatMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_float_mul", kernelInfo, comm);
  ogs::gatherManyKernel_floatMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_float_min", kernelInfo, comm);
  ogs::gatherManyKernel_floatMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_float_max", kernelInfo, comm);

  ogs::gatherManyKernel_doubleAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_double_add", kernelInfo, comm);
  ogs::gatherManyKernel_doubleMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_double_mul", kernelInfo, comm);
  ogs::gatherManyKernel_doubleMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_double_min", kernelInfo, comm);
  ogs::gatherManyKernel_doubleMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_double_max", kernelInfo, comm);

  ogs::gatherManyKernel_intAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_int_add", kernelInfo, comm);
  ogs::gatherManyKernel_intMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_int_mul", kernelInfo, comm);
  ogs::gatherManyKernel_intMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_int_min", kernelInfo, comm);
  ogs::gatherManyKernel_intMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_int_max", kernelInfo, comm);

  ogs::gatherManyKernel_longAdd = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_longlongint_add", kernelInfo, comm);
  ogs::gatherManyKernel_longMul = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_longlongint_mul", kernelInfo, comm);
  ogs::gatherManyKernel_longMin = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_longlongint_min", kernelInfo, comm);
  ogs::gatherManyKernel_longMax = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "gatherMany_longlongint_max", kernelInfo, comm);



  ogs::scatterKernel_float  = buildKernel(device, DOGS "/okl/gatherScatter.okl", "scatter_float", kernelInfo, comm);
  ogs::scatterKernel_double = buildKernel(device, DOGS "/okl/gatherScatter.okl", "scatter_double", kernelInfo, comm);
  ogs::scatterKernel_int    = buildKernel(device, DOGS "/okl/gatherScatter.okl", "scatter_int", kernelInfo, comm);
  ogs::scatterKernel_long   = buildKernel(device, DOGS "/okl/gatherScatter.okl", "scatter_long", kernelInfo, comm);

  ogs::scatterVecKernel_float  = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "scatterVec_float", kernelInfo, comm);
  ogs::scatterVecKernel_double = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "scatterVec_double", kernelInfo, comm);
  ogs::scatterVecKernel_int    = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "scatterVec_int", kernelInfo, comm);
  ogs::scatterVecKernel_long   = buildKernel(device, DOGS "/okl/gatherScatterVec.okl", "scatterVec_long", kernelInfo, comm);

  ogs::scatterManyKernel_float  = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "scatterMany_float", kernelInfo, comm);
  ogs::scatterManyKernel_double = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "scatterMany_double", kernelInfo, comm);
  ogs::scatterManyKernel_int    = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "scatterMany_int", kernelInfo, comm);
  ogs::scatterManyKernel_long   = buildKernel(device, DOGS "/okl/gatherScatterMany.okl", "scatterMany_long", kernelInfo, comm);


  if(rank==0) printf("done.\n");
}

void ogs::freeKernels() {

  ogs::gatherScatterKernel_floatAdd.free();
  ogs::gatherScatterKernel_floatMul.free();
  ogs::gatherScatterKernel_floatMin.free();
  ogs::gatherScatterKernel_floatMax.free();
  ogs::gatherScatterKernel_doubleAdd.free();
  ogs::gatherScatterKernel_doubleMul.free();
  ogs::gatherScatterKernel_doubleMin.free();
  ogs::gatherScatterKernel_doubleMax.free();
  ogs::gatherScatterKernel_intAdd.free();
  ogs::gatherScatterKernel_intMul.free();
  ogs::gatherScatterKernel_intMin.free();
  ogs::gatherScatterKernel_intMax.free();
  ogs::gatherScatterKernel_longAdd.free();
  ogs::gatherScatterKernel_longMul.free();
  ogs::gatherScatterKernel_longMin.free();
  ogs::gatherScatterKernel_longMax.free();
  ogs::gatherScatterVecKernel_floatAdd.free();
  ogs::gatherScatterVecKernel_floatMul.free();
  ogs::gatherScatterVecKernel_floatMin.free();
  ogs::gatherScatterVecKernel_floatMax.free();
  ogs::gatherScatterVecKernel_doubleAdd.free();
  ogs::gatherScatterVecKernel_doubleMul.free();
  ogs::gatherScatterVecKernel_doubleMin.free();
  ogs::gatherScatterVecKernel_doubleMax.free();
  ogs::gatherScatterVecKernel_intAdd.free();
  ogs::gatherScatterVecKernel_intMul.free();
  ogs::gatherScatterVecKernel_intMin.free();
  ogs::gatherScatterVecKernel_intMax.free();
  ogs::gatherScatterVecKernel_longAdd.free();
  ogs::gatherScatterVecKernel_longMul.free();
  ogs::gatherScatterVecKernel_longMin.free();
  ogs::gatherScatterVecKernel_longMax.free();
  ogs::gatherScatterManyKernel_floatAdd.free();
  ogs::gatherScatterManyKernel_floatMul.free();
  ogs::gatherScatterManyKernel_floatMin.free();
  ogs::gatherScatterManyKernel_floatMax.free();
  ogs::gatherScatterManyKernel_doubleAdd.free();
  ogs::gatherScatterManyKernel_doubleMul.free();
  ogs::gatherScatterManyKernel_doubleMin.free();
  ogs::gatherScatterManyKernel_doubleMax.free();
  ogs::gatherScatterManyKernel_intAdd.free();
  ogs::gatherScatterManyKernel_intMul.free();
  ogs::gatherScatterManyKernel_intMin.free();
  ogs::gatherScatterManyKernel_intMax.free();
  ogs::gatherScatterManyKernel_longAdd.free();
  ogs::gatherScatterManyKernel_longMul.free();
  ogs::gatherScatterManyKernel_longMin.free();
  ogs::gatherScatterManyKernel_longMax.free();

  ogs::gatherKernel_floatAdd.free();
  ogs::gatherKernel_floatMul.free();
  ogs::gatherKernel_floatMin.free();
  ogs::gatherKernel_floatMax.free();
  ogs::gatherKernel_doubleAdd.free();
  ogs::gatherKernel_doubleMul.free();
  ogs::gatherKernel_doubleMin.free();
  ogs::gatherKernel_doubleMax.free();
  ogs::gatherKernel_intAdd.free();
  ogs::gatherKernel_intMul.free();
  ogs::gatherKernel_intMin.free();
  ogs::gatherKernel_intMax.free();
  ogs::gatherKernel_longAdd.free();
  ogs::gatherKernel_longMul.free();
  ogs::gatherKernel_longMin.free();
  ogs::gatherKernel_longMax.free();
  ogs::gatherVecKernel_floatAdd.free();
  ogs::gatherVecKernel_floatMul.free();
  ogs::gatherVecKernel_floatMin.free();
  ogs::gatherVecKernel_floatMax.free();
  ogs::gatherVecKernel_doubleAdd.free();
  ogs::gatherVecKernel_doubleMul.free();
  ogs::gatherVecKernel_doubleMin.free();
  ogs::gatherVecKernel_doubleMax.free();
  ogs::gatherVecKernel_intAdd.free();
  ogs::gatherVecKernel_intMul.free();
  ogs::gatherVecKernel_intMin.free();
  ogs::gatherVecKernel_intMax.free();
  ogs::gatherVecKernel_longAdd.free();
  ogs::gatherVecKernel_longMul.free();
  ogs::gatherVecKernel_longMin.free();
  ogs::gatherVecKernel_longMax.free();
  ogs::gatherManyKernel_floatAdd.free();
  ogs::gatherManyKernel_floatMul.free();
  ogs::gatherManyKernel_floatMin.free();
  ogs::gatherManyKernel_floatMax.free();
  ogs::gatherManyKernel_doubleAdd.free();
  ogs::gatherManyKernel_doubleMul.free();
  ogs::gatherManyKernel_doubleMin.free();
  ogs::gatherManyKernel_doubleMax.free();
  ogs::gatherManyKernel_intAdd.free();
  ogs::gatherManyKernel_intMul.free();
  ogs::gatherManyKernel_intMin.free();
  ogs::gatherManyKernel_intMax.free();
  ogs::gatherManyKernel_longAdd.free();
  ogs::gatherManyKernel_longMul.free();
  ogs::gatherManyKernel_longMin.free();
  ogs::gatherManyKernel_longMax.free();

  ogs::scatterKernel_float.free();
  ogs::scatterKernel_double.free();
  ogs::scatterKernel_int.free();
  ogs::scatterKernel_long.free();
  ogs::scatterVecKernel_float.free();
  ogs::scatterVecKernel_double.free();
  ogs::scatterVecKernel_int.free();
  ogs::scatterVecKernel_long.free();
  ogs::scatterManyKernel_float.free();
  ogs::scatterManyKernel_double.free();
  ogs::scatterManyKernel_int.free();
  ogs::scatterManyKernel_long.free();

  ogs::o_haloBuf.free();
  ogs::haloBuf = NULL;
}

