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

#include "platform.hpp"

namespace libp {

//initialize occa::properties with common props
void platform_t::DeviceProperties(){

  properties_t& Props = props();

  Props["defines"].asObject();
  Props["includes"].asArray();
  Props["header"].asArray();
  Props["flags"].asObject();

  Props["device"].asObject();
  Props["kernel"].asObject();
  Props["memory"].asObject();

  if(sizeof(dfloat)==4){
    Props["defines/" "dfloat"]="float";
    Props["defines/" "dfloat2"]="float2";
    Props["defines/" "dfloat4"]="float4";
    Props["defines/" "dfloat8"]="float8";
  }
  if(sizeof(dfloat)==8){
    Props["defines/" "dfloat"]="double";
    Props["defines/" "dfloat2"]="double2";
    Props["defines/" "dfloat4"]="double4";
    Props["defines/" "dfloat8"]="double8";
  }

  if(sizeof(pfloat)==4){
    Props["defines/" "pfloat"]="float";
    Props["defines/" "pfloat2"]="float2";
    Props["defines/" "pfloat4"]="float4";
    Props["defines/" "pfloat8"]="float8";
  }
  if(sizeof(pfloat)==8){
    Props["defines/" "pfloat"]="double";
    Props["defines/" "pfloat2"]="double2";
    Props["defines/" "pfloat4"]="double4";
    Props["defines/" "pfloat8"]="double8";
  }

  
  if(sizeof(dlong)==4){
    Props["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
    Props["defines/" "dlong"]="long long int";
  }

  if(device.mode()=="Serial") {
    Props["compiler_flags"] += "-O3 ";
    Props["compiler_flags"] += "-g "; //debugging
    Props["defines/OCCA_USE_SERIAL"] = 1;
  }

  if(device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    Props["compiler_flags"] += "--ftz=true ";
    Props["compiler_flags"] += "--prec-div=false ";
    Props["compiler_flags"] += "--prec-sqrt=false ";
    Props["compiler_flags"] += "--use_fast_math ";
    Props["compiler_flags"] += "--fmad=true "; // compiler option for cuda
    Props["compiler_flags"] += "-Xptxas -dlcm=ca";
    Props["defines/OCCA_USE_CUDA"] = 1;
  }

  if(device.mode()=="OpenCL"){ // add backend compiler optimization for OPENCL
    Props["compiler_flags"] += " -cl-std=CL2.0 ";
    Props["compiler_flags"] += " -cl-strict-aliasing ";
    Props["compiler_flags"] += " -cl-mad-enable ";
    Props["compiler_flags"] += " -cl-no-signed-zeros ";
    Props["compiler_flags"] += " -cl-unsafe-math-optimizations ";
    Props["compiler_flags"] += " -cl-fast-relaxed-math ";
    Props["defines/OCCA_USE_OPENCL"] = 1;
  }

  if(device.mode()=="HIP"){ // add backend compiler optimization for HIP
    Props["compiler_flags"] += " -O3 ";
    Props["compiler_flags"] += " -ffp-contract=fast ";
    Props["compiler_flags"] += " -funsafe-math-optimizations ";
    Props["compiler_flags"] += " -ffast-math ";
    Props["defines/OCCA_USE_HIP"] = 1;
  }
}

} //namespace libp
