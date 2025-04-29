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

#include "stress.hpp"
#include "stressPrecon.hpp"

void stress_t::Setup(platform_t& _platform, mesh_t& _mesh,
                       settings_t& _settings, dfloat _lambda,
                       const int _NBCTypes, const memory<int> _BCType){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;
  lambda = _lambda;

  Nfields = mesh.dim;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);
  ogs::InitializeKernels(platform, ogs::Pfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"add", "sum", "scale",
        "axpy", "zaxpy",
        "amx", "amxpy", "zamxpy",
        "adx", "adxpy", "zadxpy",
        "innerProd", "norm2", "d2p", "p2d"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  // Boundary Type translation. Just defaults.
  NBCTypes = _NBCTypes;
  BCType.malloc(NBCTypes);
  BCType.copyFrom(_BCType);

  //setup boundary flags and make mask and masked ogs
  BoundarySetup();

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DSTRESS "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DSTRESS "/data/stressBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DSTRESS "/data/stressBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  
  int NblockV = std::max(1,blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  properties_t kernelInfoDouble = kernelInfo;
  kernelInfoDouble["defines/dfloat"] = "double";
  kernelInfoDouble["defines/dfloat4"] = "double4";

  properties_t kernelInfoFloat = kernelInfo;
  kernelInfoFloat["defines/dfloat"] = "float";
  kernelInfoFloat["defines/dfloat4"] = "float4";

  // Ax kernel (assume C0)
  fileName   = oklFilePrefix + "stressAx" + suffix + oklFileSuffix;
  kernelName = "stressPartialAx" + suffix;
  
  partialAxKernel = platform.buildKernel(fileName, kernelName,
					 kernelInfoDouble);
  
  floatPartialAxKernel = platform.buildKernel(fileName, kernelName,
					      kernelInfoFloat);
  
  /* Preconditioner Setup */
  Ndofs = ogsMasked.Ngather*Nfields;
  Nhalo = gHalo.Nhalo*Nfields;

  nut.malloc(mesh.Np*mesh.Nelements);
  for(int e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.Np;++n){
      dlong id = e*mesh.Np + n;
      dfloat xn = mesh.x[id];
      dfloat yn = mesh.y[id];
      nut[id] = 1 + 0.3*cos(M_PI*xn)*cos(M_PI*yn);
    }
  }
  o_nut  = platform.malloc<dfloat>(mesh.Np*mesh.Nelements, nut);


  fileName   = oklFilePrefix + "stressBuildOperatorDiagonal" + suffix + oklFileSuffix;
  kernelName = "stressBuildOperatorDiagonal" + suffix;

  if constexpr (std::is_same_v<dfloat,double>) {
    buildOperatorDiagonalKernel = platform.buildKernel(fileName, kernelName,
						       kernelInfoDouble);
  }else{
    buildOperatorDiagonalKernel = platform.buildKernel(fileName, kernelName,
						       kernelInfoFloat);
  }

  // diagonal inverse (dfloat=>(pfloat)(1/float))
  fileName   = oklFilePrefix + "stressReciprocal" + oklFileSuffix;
  kernelName = "stressReciprocal";

  if constexpr (std::is_same_v<dfloat,double>) {
    reciprocalKernel = platform.buildKernel(fileName, kernelName,
					    kernelInfoDouble);
  }else{
    reciprocalKernel = platform.buildKernel(fileName, kernelName,
					    kernelInfoFloat);
  }
  
  // assume Jacobi
  precon.Setup<JacobiPrecon>(*this);
  
}
