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

#include "mesh.hpp"
#include "mesh/mesh2D.hpp"
#include "mesh/mesh3D.hpp"

void mesh_t::MassMatrixApply(occa::memory& o_q, occa::memory& o_Mq) {
  //compute Mq = M*q
  MassMatrixKernel(Nelements, o_ggeo, o_MM, o_q, o_Mq);
}

void meshTri2D::MassMatrixKernelSetup(int Nfields) {
  occa::properties kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorTri2D.okl",
                                          "MassMatrixOperatorTri2D",
                                          kernelInfo);
}

void meshQuad2D::MassMatrixKernelSetup(int Nfields) {
  occa::properties kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorQuad2D.okl",
                                          "MassMatrixOperatorQuad2D",
                                          kernelInfo);
}

void meshTet3D::MassMatrixKernelSetup(int Nfields) {
  occa::properties kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorTet3D.okl",
                                          "MassMatrixOperatorTet3D",
                                          kernelInfo);
}

void meshHex3D::MassMatrixKernelSetup(int Nfields) {
  occa::properties kernelInfo = props; //copy base occa properties
  kernelInfo["defines/" "p_Nfields"]= Nfields;

  MassMatrixKernel = platform.buildKernel(MESH_DIR "/okl/MassMatrixOperatorHex3D.okl",
                                          "MassMatrixOperatorHex3D",
                                          kernelInfo);
}

void meshTri3D::MassMatrixKernelSetup(int Nfields) {
  LIBP_ABORT("MassMatrixOperatorTri3D not implemented yet.")
}

void meshQuad3D::MassMatrixKernelSetup(int Nfields) {
  LIBP_ABORT("MassMatrixOperatorQuad3D not implemented yet.")
}