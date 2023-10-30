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

// using namespace libp;
namespace libp {

void stab_t::stabSetupFilter(){

  int Sc=0; int Nc = 1; 

  settings.getSetting("FILTER CUTOFF", Nc);
  settings.getSetting("FILTER ORDER",  Sc);

  memory<dfloat> filterMT; 
  // Filter Settings 
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      FilterMatrixTri2D(mesh.N, Nc, Sc, filterM); 
      filterMT.malloc(filterM.length()); 
      linAlg_t::matrixTranspose(mesh.Np, mesh.Np, filterM, mesh.Np, filterMT, mesh.Np);
      break; 
    case Mesh::QUADRILATERALS:
      FilterMatrix1D(mesh.N, Nc, Sc, filterM); 
      filterMT.malloc(filterM.length()); 
      linAlg_t::matrixTranspose(mesh.Nq, mesh.Nq, filterM, mesh.Nq, filterMT, mesh.Nq);
      break; 
    case Mesh::TETRAHEDRA:
      FilterMatrixTet3D(mesh.N, Nc, Sc, filterM); 
      filterMT.malloc(filterM.length()); 
      linAlg_t::matrixTranspose(mesh.Np, mesh.Np, filterM, mesh.Np, filterMT, mesh.Np);
      break; 
    case Mesh::HEXAHEDRA:
      FilterMatrix1D(mesh.N, Nc, Sc, filterM); 
      filterMT.malloc(filterM.length()); 
      linAlg_t::matrixTranspose(mesh.Nq, mesh.Nq, filterM, mesh.Nq, filterMT, mesh.Nq);
      break; 
  }
  // Allocate Filter Matrix on Device 
  o_filterM = platform.malloc<dfloat>(filterMT); 
  
  // int blockMax = 256;
  // if (platform.device.mode() == "CUDA") blockMax = 512;

  // int sNblockV = std::max(1, blockMax/mesh.Np);
  // // int sNblockV = 1; 
  // props["defines/" "p_sNblockV"]= sNblockV;

  //  // set kernel name suffix
  //   std::string suffix;
  //   if(mesh.elementType==Mesh::TRIANGLES)
  //     suffix = "Tri2D";
  //   if(mesh.elementType==Mesh::QUADRILATERALS)
  //     suffix = "Quad2D";
  //   if(mesh.elementType==Mesh::TETRAHEDRA)
  //     suffix = "Tet3D";
  //   if(mesh.elementType==Mesh::HEXAHEDRA)
  //     suffix = "Hex3D";

  // std::string oklFilePrefix = STAB_DIR "/okl/";
  // std::string oklFileSuffix = ".okl";

  // std::string fileName, kernelName;
  // fileName      = oklFilePrefix + "filter" + suffix + oklFileSuffix;
  // kernelName    = "filter" + suffix;
  // filterKernel  = platform.buildKernel(fileName, kernelName, props);

}


// void stab_t::StabilizerApplyHJSFilter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

// // Detect Element ->o_elist!!!!
// detectorApply(o_Q, o_RHS, T); 

// // Filter solution @ o_eList
// filterKernel(mesh.Nelements, 
//              o_eList, 
//              o_filterM, 
//              o_Q); 
// }



void stab_t::FilterMatrixTri2D(int _N, int _Nc, int _sp, memory<dfloat>& _filterMatrix){
  
  // if(_Nc==0){ LIBP_FORCE_ABORT("Filter can not alter constant mode");}

  dfloat alpha = -log(std::numeric_limits<dfloat>::epsilon() );
  
  const int _Np = (_N+1)*(_N+2)/2; 
  const int _Nmodes1D = (_N+1);

  // Compute Filter Operator
  memory<dfloat> filterDiag(_Np*_Np, 0.0); 
  for(int i=0; i<_Np; i++){ filterDiag[i*_Np + i] = 1.0; }

  int sk = 0; 
  for(int i=0; i<_Nmodes1D; i++){
    for(int j=0; j<(_Nmodes1D-i); j++){
      if((i+j)>=_Nc){
         // filterDiag[sk*_Np + sk] = 0.0; 
         filterDiag[sk*_Np + sk] = exp(-alpha * pow(dfloat(i+j-_Nc)/dfloat(_N-_Nc), dfloat(_sp)));
      }
      sk++; 
    }
  }

  // Compute Filter Matrix:  V* filterdiag * invV = V*(filterdiag/V); 
  memory<dfloat> tmp(_Np*_Np), tmp_r(_Np), tmp_s(_Np); 
  mesh.NodesTri2D(_N, tmp_r, tmp_s);

  memory<dfloat> V; 
  mesh.VandermondeTri2D(_N, tmp_r, tmp_s, V);  
  
  _filterMatrix.malloc(_Np*_Np);
  linAlg_t::matrixRightSolve(_Np, _Np, filterDiag, _Np, _Np, V, tmp);
  for(int n=0; n<_Np; n++){
    for(int m=0; m<_Np; m++){
      dfloat sum = 0; 
      for(int i=0; i<_Np; i++){
        sum += V[n*_Np + i]*tmp[i*_Np + m];
      }
      _filterMatrix[n*_Np + m] = sum; 
    }
  } 

}


void stab_t::FilterMatrixTet3D(int _N, int _Nc, int _sp, memory<dfloat>& _filterMatrix){
  
  dfloat alpha = -log(std::numeric_limits<dfloat>::epsilon() );

  const int _Np = (_N+1)*(_N+2)*(_N+3)/6; 
  const int _Nmodes1D = (_N+1);

  // Compute Filter Operator
  memory<dfloat> filterDiag(_Np*_Np, 0.0); 
  for(int n=0; n<_Np; n++){ filterDiag[n*_Np + n] = 1.0; }

  int sk = 0; 
  for(int i=0; i<_Nmodes1D; i++){
    for(int j=0; j<(_Nmodes1D-i); j++){
      for(int k=0; k<(_Nmodes1D-i-j); k++){
        if((i+j+k)>=_Nc){
         // filterDiag[sk*_Np + sk] = 0.0; 
         filterDiag[sk*_Np + sk] = exp(-alpha * pow(dfloat(i+j+k -_Nc)/dfloat(_N-_Nc), dfloat(_sp)));
        }
      sk++; 
     }
   }
 }

  // Compute Filter Matrix:  V* filterdiag * invV = V*(filterdiag/V); 
  memory<dfloat> tmp(_Np*_Np), tmp_r(_Np), tmp_s(_Np), tmp_t(_Np); 
  mesh.NodesTet3D(_N, tmp_r, tmp_s, tmp_t);

  memory<dfloat> V; 
  mesh.VandermondeTet3D(_N, tmp_r, tmp_s, tmp_t, V);  
  
  _filterMatrix.malloc(_Np*_Np);
  linAlg_t::matrixRightSolve(_Np, _Np, filterDiag, _Np, _Np, V, tmp);
  for(int n=0; n<_Np; n++){
    for(int m=0; m<_Np; m++){
      dfloat sum = 0; 
      for(int i=0; i<_Np; i++){
        sum += V[n*_Np + i]*tmp[i*_Np + m];
      }
      _filterMatrix[n*_Np + m] = sum; 
    }
  } 

}


void stab_t::FilterMatrix1D(int _N, int _Nc, int _sp, memory<dfloat>& _filterMatrix){
  
  // if(_Nc==0){ LIBP_FORCE_ABORT("Filter can not alter constant mode");}

  dfloat alpha = -log(std::numeric_limits<dfloat>::epsilon() );

  const int _Np = (_N+1); 

  // Compute Filter Operator
  memory<dfloat> filterDiag(_Np*_Np, 0.0); 
  for(int i=0; i<_Np; i++){ filterDiag[i*_Np + i] = 1.0; }

  int sk = 0; 
  for(int i=0; i<_Np; i++){
      if(i>=_Nc){
         // filterDiag[sk*_Np + sk] = 0.0; 
         filterDiag[sk*_Np + sk] = exp(-alpha * pow(dfloat(i - _Nc)/ dfloat(_N -_Nc), dfloat(_sp)));
      }
      sk++; 
    }

  // Compute Filter Matrix:  V* filterdiag * invV = V*(filterdiag/V); 
  memory<dfloat> tmp(_Np*_Np), tmp_r(_Np), tmp_w(_Np); 
  mesh.JacobiGLL(_N, tmp_r, tmp_w);;

  memory<dfloat> V; 
  mesh.Vandermonde1D(_N, tmp_r, V);  
  
  _filterMatrix.malloc(_Np*_Np);
  linAlg_t::matrixRightSolve(_Np, _Np, filterDiag, _Np, _Np, V, tmp);
  for(int n=0; n<_Np; n++){
    for(int m=0; m<_Np; m++){
      dfloat sum = 0; 
      for(int i=0; i<_Np; i++){
        sum += V[n*_Np + i]*tmp[i*_Np + m];
      }
      _filterMatrix[n*_Np + m] = sum; 
    }
  } 

}


} //namespace libp
