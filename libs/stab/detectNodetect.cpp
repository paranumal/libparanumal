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

#include "stab.hpp"

namespace libp {

void stab_t::detectSetupNodetect(){
  // Initialize Required Memory
  elementList.malloc((mesh.Nelements+mesh.totalHaloPairs)*Ndfields); 
  o_elementList = platform.malloc<dlong>(elementList); 

  // Create a field for detector
  qdetector.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*Ndfields); 
  o_qdetector = platform.malloc<dfloat>(qdetector); 

  // Compute Art. Diff. Activation function
  if(type==Stab::ARTDIFF){ // Vertex based to make it continuous
    viscosityActivation.malloc((mesh.Nelements+mesh.totalHaloPairs)*Ndfields, 0.0); 
    o_viscosityActivation = platform.malloc<dfloat>(viscosityActivation); 
  }
   props["defines/" "p_blockSize"]= 256;
   props["defines/" "s_Ndfields"]= Ndfields;
   props["defines/" "s_Nsfields"]= Nsfields;
  // Needed for Subcell and Limiter Only
  if(type==Stab::SUBCELL || type==Stab::LIMITER){
     props["defines/" "s_DGDG_TYPE"] = int(0); 
     props["defines/" "s_FVFV_TYPE"] = int(1); 
     props["defines/" "s_DGFV_TYPE"] = int(2); 
  }  
}


void stab_t::detectApplyNodetect(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
platform.linAlg().set((mesh.Nelements+mesh.totalHaloPairs)*Ndfields, 0, o_elementList);
platform.linAlg().set((mesh.Nelements+mesh.totalHaloPairs)*Ndfields, 0.0, o_viscosityActivation);
}

} //namespace libp
