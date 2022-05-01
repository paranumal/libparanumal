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

#include "elliptic.hpp"
#include "ellipticPrecon.hpp"

void elliptic_t::Setup(platform_t& _platform, mesh_t& _mesh,
                       settings_t& _settings, dfloat _lambda,
                       const int _NBCTypes, const memory<int> _BCType){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;
  lambda = _lambda;

  Nfields = 1;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  disc_ipdg = settings.compareSetting("DISCRETIZATION","IPDG");
  disc_c0   = settings.compareSetting("DISCRETIZATION","CONTINUOUS");

  //setup linear algebra module
  platform.linAlg().InitKernels({"add", "sum", "scale",
                                "axpy", "zaxpy",
                                "amx", "amxpy", "zamxpy",
                                "adx", "adxpy", "zadxpy",
                                "innerProd", "norm2"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // Boundary Type translation. Just defaults.
  NBCTypes = _NBCTypes;
  BCType.malloc(NBCTypes);
  BCType.copyFrom(_BCType);

  //setup boundary flags and make mask and masked ogs
  BoundarySetup();

  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    //tau (penalty term in IPDG)
    if (mesh.elementType==Mesh::TRIANGLES ||
        mesh.elementType==Mesh::QUADRILATERALS){
      tau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      if(mesh.dim==3)
        tau *= 1.5;
    } else
      tau = 2.0*(mesh.N+1)*(mesh.N+3);

    //buffer for gradient
    dlong Ntotal = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
    grad.malloc(Ntotal*4);
    o_grad = platform.malloc<dfloat>(grad);
  } else {
    tau = 0.0;

    //buffer for local Ax
    dlong Ntotal = mesh.Np*mesh.Nelements;
    o_AqL = platform.malloc<dfloat>(Ntotal);
  }

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES){
    if(mesh.dim==2)
      suffix = "Tri2D";
    else
      suffix = "Tri3D";
  } else if(mesh.elementType==Mesh::QUADRILATERALS){
    if(mesh.dim==2)
      suffix = "Quad2D";
    else
      suffix = "Quad3D";
  } else if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  else if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";

  std::string oklFilePrefix = DELLIPTIC "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  // Ax kernel
  if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    fileName   = oklFilePrefix + "ellipticAx" + suffix + oklFileSuffix;
    if(mesh.elementType==Mesh::HEXAHEDRA){
      if(mesh.settings.compareSetting("ELEMENT MAP", "TRILINEAR"))
        kernelName = "ellipticPartialAxTrilinear" + suffix;
      else
        kernelName = "ellipticPartialAx" + suffix;
    } else{
      kernelName = "ellipticPartialAx" + suffix;
    }

    partialAxKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);

  } else if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    int Nmax = std::max(mesh.Np, mesh.Nfaces*mesh.Nfp);
    kernelInfo["defines/" "p_Nmax"]= Nmax;

    fileName   = oklFilePrefix + "ellipticGradient" + suffix + oklFileSuffix;
    kernelName = "ellipticPartialGradient" + suffix;
    partialGradientKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

    fileName   = oklFilePrefix + "ellipticAxIpdg" + suffix + oklFileSuffix;
    kernelName = "ellipticPartialAxIpdg" + suffix;
    partialIpdgKernel = platform.buildKernel(fileName, kernelName,
                                              kernelInfo);
  }

  /* Preconditioner Setup */
  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    Ndofs = ogsMasked.Ngather*Nfields;
    Nhalo = gHalo.Nhalo*Nfields;
  } else {
    Ndofs = mesh.Nelements*mesh.Np*Nfields;
    Nhalo = mesh.totalHaloPairs*mesh.Np*Nfields;
  }

  if       (settings.compareSetting("PRECONDITIONER", "JACOBI"))
    precon.Setup<JacobiPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "MASSMATRIX"))
    precon.Setup<MassMatrixPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "PARALMOND"))
    precon.Setup<ParAlmondPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "MULTIGRID"))
    precon.Setup<MultiGridPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "SEMFEM"))
    precon.Setup<SEMFEMPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "OAS"))
    precon.Setup<OASPrecon>(*this);
  else if(settings.compareSetting("PRECONDITIONER", "NONE"))
    precon.Setup<IdentityPrecon>(Ndofs);
}
