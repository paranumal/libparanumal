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

#include "mesh.hpp"

namespace libp {

void mesh_t::MassMatrixApply(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_Mq) {
  //compute Mq = M*q
  MassMatrixKernel(Nelements, o_wJ, o_MM, o_q, o_Mq);
}

void mesh_t::MassMatrixKernelSetupTri2D(int Nfields) {
  properties_t kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorTri2D.okl",
                                          "MassMatrixOperatorTri2D",
                                          kernelInfo);
}

void mesh_t::MassMatrixKernelSetupQuad2D(int Nfields) {
  properties_t kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorQuad2D.okl",
                                          "MassMatrixOperatorQuad2D",
                                          kernelInfo);
}

void mesh_t::MassMatrixKernelSetupTet3D(int Nfields) {
  properties_t kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorTet3D.okl",
                                          "MassMatrixOperatorTet3D",
                                          kernelInfo);
}

void mesh_t::MassMatrixKernelSetupHex3D(int Nfields) {
  properties_t kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorHex3D.okl",
                                          "MassMatrixOperatorHex3D",
                                          kernelInfo);
}

} //namespace libp
