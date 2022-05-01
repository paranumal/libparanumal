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

#include "gradient.hpp"

void gradient_t::Setup(platform_t& _platform, mesh_t& _mesh,
                       gradientSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = mesh.comm;
  settings = _settings;

  Nfields = mesh.dim;

  dlong Nlocal = mesh.Nelements*mesh.Np;

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal);
  o_q = platform.malloc<dfloat>(q);

  gradq.malloc(Nlocal*mesh.dim);
  o_gradq = platform.malloc<dfloat>(gradq);

  //storage for M*gradq during reporting
  o_Mgradq = platform.malloc<dfloat>(gradq);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

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

  std::string oklFilePrefix = DGRADIENT "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "gradientVolume" + suffix + oklFileSuffix;
  kernelName = "gradientVolume" + suffix;

  volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "gradientInitialCondition2D" + oklFileSuffix;
    kernelName = "gradientInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "gradientInitialCondition3D" + oklFileSuffix;
    kernelName = "gradientInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

}
