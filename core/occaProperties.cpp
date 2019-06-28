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

#include <occa.hpp>
#include "types.h"
#include "utils.hpp"

//initialize occa::properties with common props
void occaDeviceProperties(occa::device &device, occa::properties& props){

  props["defines"].asObject();
  props["includes"].asArray();
  props["header"].asArray();
  props["flags"].asObject();

  props["device"].asObject();
  props["kernel"].asObject();
  props["memory"].asObject();

  if(sizeof(dfloat)==4){
    props["defines/" "dfloat"]="float";
    props["defines/" "dfloat2"]="float2";
    props["defines/" "dfloat4"]="float4";
    props["defines/" "dfloat8"]="float8";
  }
  if(sizeof(dfloat)==8){
    props["defines/" "dfloat"]="double";
    props["defines/" "dfloat2"]="double2";
    props["defines/" "dfloat4"]="double4";
    props["defines/" "dfloat8"]="double8";
  }

  if(sizeof(dlong)==4){
    props["defines/" "dlong"]="int";
  }
  if(sizeof(dlong)==8){
    props["defines/" "dlong"]="long long int";
  }

  if(device.mode()=="Serial")
    props["compiler_flags"] += "-g"; //debugging

  if(device.mode()=="CUDA"){ // add backend compiler optimization for CUDA
    props["compiler_flags"] += "--ftz=true ";
    props["compiler_flags"] += "--prec-div=false ";
    props["compiler_flags"] += "--prec-sqrt=false ";
    props["compiler_flags"] += "--use_fast_math ";
    props["compiler_flags"] += "--fmad=true "; // compiler option for cuda
    props["compiler_flags"] += "-Xptxas -dlcm=ca";
  }
}
