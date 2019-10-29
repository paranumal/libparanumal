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

#include "gradient.hpp"

gradient_t& gradient_t::Setup(mesh_t& mesh, linAlg_t& linAlg){

  gradient_t* gradient = new gradient_t(mesh, linAlg);

  settings_t& settings = gradient->settings;

  gradient->Nfields = mesh.dim;

  dlong Nlocal = mesh.Nelements*mesh.Np;

  //setup linear algebra module
  gradient->linAlg.InitKernels({"innerProd"}, mesh.comm);

  // compute samples of q at interpolation nodes
  gradient->q = (dfloat*) calloc(Nlocal, sizeof(dfloat));
  gradient->o_q = mesh.device.malloc(Nlocal*sizeof(dfloat), gradient->q);

  gradient->gradq = (dfloat*) calloc(Nlocal*mesh.dim, sizeof(dfloat));
  gradient->o_gradq = mesh.device.malloc(Nlocal*mesh.dim*sizeof(dfloat), gradient->gradq);

  //storage for M*gradq during reporting
  gradient->o_Mgradq = mesh.device.malloc(Nlocal*mesh.dim*sizeof(dfloat), gradient->gradq);

  // OCCA build stuff
  occa::properties kernelInfo = gradient->props; //copy base occa properties

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= gradient->Nfields;

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  // kernels from volume file
  sprintf(fileName, DGRADIENT "/okl/gradientVolume%s.okl", suffix);
  sprintf(kernelName, "gradientVolume%s", suffix);

  gradient->volumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  // mass matrix operator
  sprintf(fileName, LIBP_DIR "/core/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  gradient->MassMatrixKernel = buildKernel(mesh.device, fileName, kernelName,
                                            kernelInfo, mesh.comm);


  if (mesh.dim==2) {
    sprintf(fileName, DGRADIENT "/okl/gradientInitialCondition2D.okl");
    sprintf(kernelName, "gradientInitialCondition2D");
  } else {
    sprintf(fileName, DGRADIENT "/okl/gradientInitialCondition3D.okl");
    sprintf(kernelName, "gradientInitialCondition3D");
  }

  gradient->initialConditionKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);

  return *gradient;
}

gradient_t::~gradient_t() {
  volumeKernel.free();
  MassMatrixKernel.free();
  initialConditionKernel.free();
}