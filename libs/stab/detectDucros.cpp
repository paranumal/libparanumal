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

void stab_t::detectSetupDucros(){
  /*
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

  // Create required operators 
  memory<dfloat> invV, invVT;
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      ModeInfoKlocknerTri2D(mesh.N, modeMap); 
      mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, invV);
      break; 
    case Mesh::QUADRILATERALS:
      ModeInfoKlocknerQuad2D(mesh.N, modeMap); 
      mesh.VandermondeQuad2D(mesh.N, mesh.r, mesh.s, invV);
      break; 
    case Mesh::TETRAHEDRA:
      ModeInfoKlocknerTet3D(mesh.N, modeMap); 
      mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, invV);
      break; 
    case Mesh::HEXAHEDRA:
      ModeInfoKlocknerHex3D(mesh.N, modeMap); 
      mesh.VandermondeHex3D(mesh.N, mesh.r, mesh.s, mesh.t, invV);
      break; 
  }

   o_modeMap = platform.malloc<int>(modeMap); 
  
  // 
  invVT.malloc(mesh.Np*mesh.Np, 0.0); 
  linAlg_t::matrixInverse(mesh.Np, invV); 
  linAlg_t::matrixTranspose(mesh.Np, mesh.Np, invV, mesh.Np, invVT, mesh.Np);

   o_invV    = platform.malloc<dfloat>(invVT); 

  // Klockner 1D mode operations 
   LeastSquaresFitKlockner(mesh.N, leastSquares1D); 
   BaseLineDecayKlockner(mesh.N, baseLineDecay); 
   o_leastSquares1D = platform.malloc<dfloat>(leastSquares1D); 
   o_baseLineDecay  = platform.malloc<dfloat>(baseLineDecay); 

   props["defines/" "p_blockSize"]= 256;
   props["defines/" "s_Ndfields"]= Ndfields;
   props["defines/" "s_Nsfields"]= Nsfields;
   props["defines/" "s_sS0"] = mesh.N <= 3 ? 1.0:2.0;  
   props["defines/" "s_sK0"] = 1.0;  
   props["defines/" "p_Nq"]  = mesh.N+1;

   // props["defines/" "s_Nverts"]  = mesh.Nverts;  

  // Needed for Subcell and Limiter Only
  if(type==Stab::SUBCELL || type==Stab::LIMITER){
     props["defines/" "s_DGDG_TYPE"] = int(0); 
     props["defines/" "s_FVFV_TYPE"] = int(1); 
     props["defines/" "s_DGFV_TYPE"] = int(2); 
  }

 // set kernel name suffix
  std::string suffix        = mesh.elementSuffix();
  std::string oklFilePrefix = STAB_DIR "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName      = oklFilePrefix + "detect" + suffix + oklFileSuffix;

  // if artificial diffusion return activation function
  if(type==Stab::ARTDIFF){
    kernelName    = "detectKlocknerDiffusion" + suffix;
  }else{ // return coefficient
    kernelName    = "detectKlockner" + suffix;
  }

  detectKernel  = platform.buildKernel(fileName, kernelName, props);

  if(type==Stab::SUBCELL || type==Stab::LIMITER){
    fileName        = oklFilePrefix + "subcell" + oklFileSuffix;
    kernelName      = "detectFindNeigh" + suffix;
    findNeighKernel =  platform.buildKernel(fileName, kernelName, props);
  }
  
  fileName        = oklFilePrefix + "utilities" + oklFileSuffix; 

  kernelName         = "extractField";
  extractFieldKernel = platform.buildKernel(fileName, kernelName, props); 

  // fileName        = oklFilePrefix + "utilities" + oklFileSuffix; 
  // kernelName      = "copyFloat";
  // copyFloatKernel = platform.buildKernel(fileName, kernelName, props);  

  // kernelName      = "copyInt";
  // copyIntKernel   = platform.buildKernel(fileName, kernelName, props); 

  
  // // Initialize Required Memeory: Remove LaterAK!
  // efList.malloc((mesh.Nelements+mesh.totalHaloPairs)*dNfields); 
  // o_efList = platform.malloc<dfloat>(efList); 

*/
  
}


void stab_t::detectApplyDucros(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
/*
// // Directly copy from field for HJS
// if(solver==Stab::HJS){
  o_qdetector.copyFrom(o_Q); 
// }else if(solver==Stab::CNS){
//   // Use Density field for now AK:
//   const int field_id = 0; 
//   extractFieldKernel(mesh.Nelements,
//                      field_id, 
//                      o_Q,
//                      o_qdetector); 
// }

if(type==Stab::ARTDIFF){
  // Detect elements for each fields i.e. 2
  detectKernel(mesh.Nelements, 
               mesh.o_vgeo, 
               o_modeMap, 
               mesh.o_MM, 
               o_invV, 
               o_leastSquares1D, 
               o_baseLineDecay, 
               o_qdetector, 
               o_viscosityActivation,
               o_elementList); 
  
  }else if(type==Stab::SUBCELL || 
           type==Stab::LIMITER || 
           type==Stab::FILTER){

  detectKernel(mesh.Nelements, 
               mesh.o_vgeo, 
               o_modeMap, 
               mesh.o_MM, 
               o_invV, 
               o_leastSquares1D, 
               o_baseLineDecay, 
               o_qdetector, 
               o_elementList); 

  // make sure that this works for int!!!!!! AK.
  mesh.halo.Exchange(o_elementList, Ndfields); 

  // Correct Neighbor Info
  findNeighKernel(mesh.Nelements, 
                  mesh.o_vmapP, 
                  o_elementList); 

}
*/
}




} //namespace libp
