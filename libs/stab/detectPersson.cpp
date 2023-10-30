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

void stab_t::detectSetupPersson(){

  // Initialize Required Memeory
  elementList.malloc((mesh.Nelements+mesh.totalHaloPairs)*Ndfields); 
  o_elementList = platform.malloc<dlong>(elementList); 
/*
  // Initialize Required Memeory: Remove LaterAK!
  efList.malloc((mesh.Nelements+mesh.totalHaloPairs)*dNfields); 
  o_efList = platform.malloc<dfloat>(efList); 

  qd.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*dNfields); 
  o_qd = platform.malloc<dfloat>(qd); 


  // Project to modal space spans N-1
  projectNm1.malloc(mesh.Np*mesh.Np);

  // Compute Art. Diff. Activation function as well
  if(type==Stab::ARTDIFF){ // Vertex based to make it continuous
    viscActivation.malloc((mesh.Nelements+mesh.totalHaloPairs)*dNfields, 0.0); 
    o_viscActivation = platform.malloc<dfloat>(viscActivation); 
  }


  memory<dfloat> V, invV, tmp, truncModes;
  // Compute 2D to 1D mode map 
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      ModeInfoPerssonTri2D(mesh.N, truncModes); 
      mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, V);
      break; 
    case Mesh::QUADRILATERALS:
      ModeInfoPerssonQuad2D(mesh.N, truncModes); 
      mesh.VandermondeQuad2D(mesh.N, mesh.r, mesh.s, V);
      break; 
    case Mesh::TETRAHEDRA:
      ModeInfoPerssonTet3D(mesh.N, truncModes); 
      mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
      break; 
    case Mesh::HEXAHEDRA:
      ModeInfoPerssonHex3D(mesh.N, truncModes); 
      mesh.VandermondeHex3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
      break; 
  }

  // Copy Vandermonde and invert
  invV.malloc(mesh.Np*mesh.Np); 
  invV.copyFrom(V, mesh.Np*mesh.Np); // = V.clone(); 
  linAlg_t::matrixInverse(mesh.Np, invV); 

  // Cut-off Highest Modes
  tmp.malloc(mesh.Np*mesh.Np); 

  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Np; j++){
      tmp[i*mesh.Np + j] = truncModes[i]*invV[i*mesh.Np + j]; 
    }
  }

  // Transponse of Projection Operator
  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Np; j++){
      dfloat sum = 0; 
      for(int m=0; m<mesh.Np; m++){
          sum += V[i*mesh.Np +m]*tmp[m*mesh.Np + j]; 
      }
      projectNm1[j*mesh.Np + i] = sum; 
    }
  }

   o_projectNm1    = platform.malloc<dfloat>(projectNm1); 

   props["defines/" "p_dNfields"]= dNfields;
   props["defines/" "p_sNfields"]= sNfields;
   props["defines/" "p_sNverts"] = mesh.Nverts;  
   props["defines/" "p_Nq"]= mesh.N+1;

    // Needed for Subcell Only
    if(type==Stab::SUBCELL || type==Stab::LIMITER){
       props["defines/" "s_DGDG_TYPE"] = int(0); 
       props["defines/" "s_FVFV_TYPE"] = int(1); 
       props["defines/" "s_DGFV_TYPE"] = int(2); 
    }

   // set kernel name suffix
    std::string suffix;
    if(mesh.elementType==Mesh::TRIANGLES)
      suffix = "Tri2D";
    if(mesh.elementType==Mesh::QUADRILATERALS)
      suffix = "Quad2D";
    if(mesh.elementType==Mesh::TETRAHEDRA)
      suffix = "Tet3D";
    if(mesh.elementType==Mesh::HEXAHEDRA)
      suffix = "Hex3D";

  std::string oklFilePrefix = STAB_DIR "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  fileName      = oklFilePrefix + "detect" + suffix + oklFileSuffix;
  if(type==Stab::ARTDIFF){
    kernelName    = "detectPerssonDiffusion" + suffix;
  }else{
    kernelName    = "detectPersson" + suffix;
  }
  detectKernel  = platform.buildKernel(fileName, kernelName, props);

  if(type==Stab::SUBCELL || type==Stab::LIMITER){
    fileName      = oklFilePrefix + "subcell" + oklFileSuffix;
    kernelName    = "detectFindNeigh" + suffix;
    findNeighKernel =  platform.buildKernel(fileName, kernelName, props);
  }

  fileName        = oklFilePrefix + "utilities" + oklFileSuffix; 
  kernelName      = "copyFloat";
  copyFloatKernel = platform.buildKernel(fileName, kernelName, props);  

  kernelName      = "copyInt";
  copyIntKernel   = platform.buildKernel(fileName, kernelName, props); 

  kernelName         = "extractField";
  extractFieldKernel = platform.buildKernel(fileName, kernelName, props); 
  */
}


void stab_t::detectApplyPersson(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){
/*
// Directly copy from field for HJS
if(solver==Stab::HJS){
  o_qd.copyFrom(o_Q); 
}else if(solver==Stab::CNS){
  // Use Density field for now AK:
  const int field_id = 0; 
  extractFieldKernel(mesh.Nelements,
                     field_id, 
                     o_Q,
                     o_qd); 
}

if(type==Stab::ARTDIFF){
  // Detect elements for each fields i.e. 2
  detectKernel(mesh.Nelements, 
               mesh.o_vgeo, 
               mesh.o_MM, 
               o_projectNm1, 
               o_qd,
               o_viscActivation, 
               o_eList);


   }else if(type==Stab::SUBCELL || type==Stab::LIMITER){
   // Detect elements for each fields i.e. 2
  detectKernel(mesh.Nelements, 
               mesh.o_vgeo, 
               mesh.o_MM, 
               o_projectNm1, 
               o_qd,
               o_efList); 

  // Communicate as we need neighbor info
  mesh.halo.Exchange(o_efList, dNfields); 

  // Correct Neighbor Info
  findNeighKernel(mesh.Nelements, 
                  mesh.o_vmapP, 
                  o_efList); 

  // Return to int storage
  copyIntKernel((mesh.Nelements+mesh.totalHaloPairs)*dNfields, 
              o_efList, 
              o_eList); 
}else{
    // Detect elements for each fields i.e. 2
    detectKernel(mesh.Nelements, 
                 mesh.o_vgeo, 
                 mesh.o_MM, 
                 o_projectNm1, 
                 o_qd,
                 o_eList); 

 }


 // mesh.halo.Exchange(o_eList, dNfields); 

 // printf("%d\t%d\n", mesh.NhaloElements, mesh.totalHaloPairs);
*/
}




void stab_t::ModeInfoPerssonTri2D(int _N, memory<dfloat>& _truncModes){
  const int _Np = (_N+1)*(_N+2)/2; 
  const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        if(i+j < (_Nmodes1D-1)){
          _truncModes[sk++] = 1.;
        }else{ 
          _truncModes[sk++] = 0.;
        }
      }
    }
}


void stab_t::ModeInfoPerssonQuad2D(int _N, memory<dfloat>& _truncModes){
  const int _Np = (_N+1)*(_N+1); 
  const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D; j++){
        if(std::max(i,j) < (_Nmodes1D-1)){
          _truncModes[sk++] = 1.;
        }else{ 
          _truncModes[sk++] = 0.;
        }
      }
    }
}



void stab_t::ModeInfoPerssonTet3D(int _N, memory<dfloat>& _truncModes){
  const int _Np = (_N+1)*(_N+2)*(_N+3)/6; 
  const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        for(int k=0; k<_Nmodes1D-i-j; k++){
          if((i+j+k)< (_Nmodes1D-1)){
           _truncModes[sk++] = 1.;
            }else{ 
          _truncModes[sk++] = 0.;
          }
        }
      }
    }
}

void stab_t::ModeInfoPerssonHex3D(int _N, memory<dfloat>& _truncModes){
 const int _Np = (_N+1)*(_N+1)*(_N+1); 
 const int _Nmodes1D = (_N+1); 
  _truncModes.malloc(_Np); 

  int sk = 0; 
    for(int j=0; j<_Nmodes1D; j++){
      for(int i=0; i<_Nmodes1D; i++){
        for(int k=0; k<_Nmodes1D; k++){
         if(std::max(std::max(i,j),k) < (_Nmodes1D-1)){
            _truncModes[sk++] = 1.;
            }else{ 
            _truncModes[sk++] = 0.;
          }
         }
        }
      }
  }



} //namespace libp
