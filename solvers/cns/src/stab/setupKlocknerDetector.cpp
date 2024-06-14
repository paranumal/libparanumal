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

#include "cns.hpp"

void ModeInfoKlocknerTri2D(int _N, memory<int>& _modeMap);
void ModeInfoKlocknerQuad2D(int _N, memory<int>&_modeMap);
void ModeInfoKlocknerTet3D(int _N, memory<int>&_modeMap);
void ModeInfoKlocknerHex3D(int _N, memory<int>&_modeMap);
void LeastSquaresFit(int _N, memory<dfloat>& _LSF);
void BaseLineDecay(int _N, memory<dfloat>& _BLD);


void cns_t::applyKlocknerDetector(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T){

  detectKernel(mesh.Nelements, 
              mesh.o_vgeo, 
              o_modeIds, 
              mesh.o_MM, 
              o_invV, 
              o_LS1D, 
              o_BLD, 
              o_Q, 
              o_qdetect, 
              o_elmList);
}


void cns_t::setupKlocknerDetector(){

  // There is no cutout in the detector to mark elements
  // as this should be function of stabilizer 
  // Initialize Required Memory
  elmList.malloc((mesh.Nelements+mesh.totalHaloPairs)*1); 
  o_elmList = platform.malloc<dlong>(elmList); 

  qdetect.malloc((mesh.Nelements+mesh.totalHaloPairs)*1); 
  o_qdetect = platform.malloc<dfloat>(qdetect); 

   // Create required operators 
  memory<int> modeIds; 
  memory<dfloat> invV, invVT;
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      ModeInfoKlocknerTri2D(mesh.N, modeIds); 
      mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, invV);
      break; 
    case Mesh::QUADRILATERALS:
      ModeInfoKlocknerQuad2D(mesh.N, modeIds); 
      mesh.VandermondeQuad2D(mesh.N, mesh.r, mesh.s, invV);
      break; 
    case Mesh::TETRAHEDRA:
      ModeInfoKlocknerTet3D(mesh.N, modeIds); 
      mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, invV);
      break; 
    case Mesh::HEXAHEDRA:
      ModeInfoKlocknerHex3D(mesh.N, modeIds); 
      mesh.VandermondeHex3D(mesh.N, mesh.r, mesh.s, mesh.t, invV);
      break; 
  }

  o_modeIds = platform.malloc<int>(modeIds); 

   // 
  invVT.malloc(mesh.Np*mesh.Np, 0.0); 
  linAlg_t::matrixInverse(mesh.Np, invV); 
  linAlg_t::matrixTranspose(mesh.Np, mesh.Np, invV, mesh.Np, invVT, mesh.Np);
  o_invV    = platform.malloc<dfloat>(invVT); 

  // Klockner 1D mode operations
  memory<dfloat> LS1D, BLD;  
  LeastSquaresFit(mesh.N, LS1D); 
  BaseLineDecay(mesh.N, BLD); 
  o_LS1D = platform.malloc<dfloat>(LS1D); 
  o_BLD  = platform.malloc<dfloat>(BLD);

  props["defines/" "s_sS0"] = mesh.N <= 3 ? 1.0:2.0;  
  // props["defines/" "s_sS0"] = 2.0; 
  props["defines/" "s_sK0"] = 1.0;  
  props["defines/" "p_Nq"]  = mesh.N+1;

   // set kernel name suffix
  std::string suffix        = mesh.elementSuffix();
  std::string oklFilePrefix = DCNS"/okl/detect/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;
  fileName      = oklFilePrefix + "cnsDetectKlockner" + oklFileSuffix;
  kernelName    = "cnsDetectKlockner" + suffix;
  detectKernel  = platform.buildKernel(fileName, kernelName, props);

}


// Utility Functions

// Create ModeMap for increasing order i.e. [[0], [1, N+1], [2 N+2 2*N+1], etc.]
// such that [zero-order, first-order, second-order]
void ModeInfoKlocknerTri2D(int _N, memory<int>& _modeMap){
  const int _Np = (_N+1)*(_N+2)/2; 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        if(i+j == id){
          _modeMap[n++] = sk; 
        }
      sk++;
          }
        }
    sk=0;
  }

}

void ModeInfoKlocknerQuad2D(int _N, memory<int>&_modeMap){
  const int _Np = (_N+1)*(_N+1); 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int j=0; j<_Nmodes1D; j++){
      for(int i=0; i<_Nmodes1D; i++){
         if(std::max(i,j) == id){
            _modeMap[n++] = sk; 
          }
        sk++;
         }
       }
    sk=0;
  }
}

void ModeInfoKlocknerTet3D(int _N, memory<int>& _modeMap){
  const int _Np = (_N+1)*(_N+2)*(_N+3)/6; 
  const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int i=0; i<_Nmodes1D; i++){
      for(int j=0; j<_Nmodes1D-i; j++){
        for(int k=0; k<_Nmodes1D-i-j; k++){
          if((i+j+k)==id){
            _modeMap[n++] = sk;
          } 
          sk++;
        }
      }
    }
    sk=0;
  }
}



void ModeInfoKlocknerHex3D(int _N, memory<int>& _modeMap){
 const int _Np = (_N+1)*(_N+1)*(_N+1); 
 const int _Nmodes1D = (_N+1); 
  _modeMap.malloc(_Np); 

  int sk = 0, n=0; 
  for(int id=0; id<_Nmodes1D;id++){
    for(int j=0; j<_Nmodes1D; j++){
      for(int i=0; i<_Nmodes1D; i++){
        for(int k=0; k<_Nmodes1D; k++){
         if(std::max(std::max(i,j),k) == id){
            _modeMap[n++] = sk; 
          }
        sk++;
         }
        }
      }
    sk=0;
  }
}



// Least squares fit for 1D modes
void LeastSquaresFit(int _N, memory<dfloat>& _LSF){
  memory<dfloat>  tmp(2*_N); 
  _LSF.malloc( _N); 

  for(int n=0; n<_N; n++){
    const dfloat logmode = log10(n+1); 
    tmp[2*n + 0] = logmode; 
    tmp[2*n + 1] = 1.0; 
  }

  linAlg_t::matrixPseudoInverse(_N, 2, tmp);

  for(int n=0; n<_N; n++){
    _LSF[n] = tmp[n];
  }
}


//  baseline decay (squared)
void BaseLineDecay(int _N, memory<dfloat>& _BLD){
  _BLD.malloc(_N+1, 0.0);   
  dfloat bsum = 0.0; 
  for(int j=1; j<_N+1; j++)
    bsum +=1.0/pow(j, 2*_N); 

  bsum = 1.0/sqrt(bsum); 

  _BLD[0] = 0.0; 
  // baseline decay (squared) 
  for(int n=1; n<_N+1; n++){
    const dfloat bdecay = bsum*1.0/(pow(n,_N));
    _BLD[n] = bdecay*bdecay;
  }
}


