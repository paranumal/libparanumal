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

#include "elliptic.hpp"
#include "ellipticPrecon.hpp"

elliptic_t& elliptic_t::Setup(mesh_t& mesh, linAlg_t& linAlg, dfloat lambda){

  elliptic_t* elliptic = new elliptic_t(mesh, linAlg, lambda);

  settings_t& settings = elliptic->settings;

  elliptic->Nfields = 1;

  elliptic->disc_ipdg = settings.compareSetting("DISCRETIZATION","IPDG");
  elliptic->disc_c0   = settings.compareSetting("DISCRETIZATION","CONTINUOUS");

  //setup linear algebra module
  elliptic->linAlg.InitKernels({"add", "sum", "scale",
                                "axpy", "zaxpy",
                                "amx", "amxpy", "zamxpy",
                                "adx", "adxpy", "zadxpy",
                                "innerProd", "weightedInnerProd",
                                "norm2", "weightedNorm2"},
                                mesh.comm);

  /*setup trace halo exchange */
  elliptic->traceHalo = mesh.HaloTraceSetup(elliptic->Nfields);

  // Boundary Type translation. Just defaults.
  int BCType[3] = {0,1,2};
  elliptic->BCType = (int*) calloc(3,sizeof(int));
  memcpy(elliptic->BCType,BCType,3*sizeof(int));

  //setup boundary flags and make mask and masked ogs
  elliptic->BoundarySetup();

  //tau (penalty term in IPDG)
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    if (mesh.elementType==TRIANGLES ||
        mesh.elementType==QUADRILATERALS){
      elliptic->tau = 2.0*(mesh.N+1)*(mesh.N+2)/2.0;
      if(mesh.dim==3)
        elliptic->tau *= 1.5;
    } else
      elliptic->tau = 2.0*(mesh.N+1)*(mesh.N+3);

    //buffer for gradient
    dlong Ntotal = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
    elliptic->grad = (dfloat*) calloc(Ntotal*4, sizeof(dfloat));
    elliptic->o_grad  = mesh.device.malloc(Ntotal*4*sizeof(dfloat), elliptic->grad);
  }

  // OCCA build stuff
  occa::properties kernelInfo = elliptic->props; //copy base occa properties

  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES){
    if(mesh.dim==2)
      suffix = strdup("Tri2D");
    else
      suffix = strdup("Tri3D");
  } else if(mesh.elementType==QUADRILATERALS){
    if(mesh.dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  } else if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  else if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  // mask
  elliptic->maskKernel = buildKernel(mesh.device, DELLIPTIC "/okl/ellipticMask.okl",
                                     "mask", kernelInfo, mesh.comm);

  //add standard boundary functions
  char *boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int NblockV = mymax(1,512/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  // Ax kernel
  if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
    if(mesh.elementType==HEXAHEDRA){
      if(settings.compareSetting("ELEMENT MAP", "TRILINEAR"))
        sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
      else
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
    } else{
      sprintf(kernelName, "ellipticPartialAx%s", suffix);
    }

    elliptic->partialAxKernel = buildKernel(mesh.device, fileName, kernelName,
                                     kernelInfo, mesh.comm);

  } else if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    int Nmax = mymax(mesh.Np, mesh.Nfaces*mesh.Nfp);
    kernelInfo["defines/" "p_Nmax"]= Nmax;

    sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
    sprintf(kernelName, "ellipticPartialGradient%s", suffix);
    elliptic->partialGradientKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);

    sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
    sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
    elliptic->partialIpdgKernel = buildKernel(mesh.device, fileName, kernelName,
                                              kernelInfo, mesh.comm);
  }

  /* Preconditioner Setup */
  if       (settings.compareSetting("PRECONDITIONER", "JACOBI"))
    elliptic->precon = new JacobiPrecon(*elliptic);
  else if(settings.compareSetting("PRECONDITIONER", "MASSMATRIX"))
    elliptic->precon = new MassMatrixPrecon(*elliptic);
  else if(settings.compareSetting("PRECONDITIONER", "FULLALMOND"))
    elliptic->precon = new ParAlmondPrecon(*elliptic);
  else if(settings.compareSetting("PRECONDITIONER", "MULTIGRID"))
    elliptic->precon = new MultiGridPrecon(*elliptic);
  else if(settings.compareSetting("PRECONDITIONER", "SEMFEM"))
    elliptic->precon = new SEMFEMPrecon(*elliptic);
  else if(settings.compareSetting("PRECONDITIONER", "OAS"))
    elliptic->precon = new OASPrecon(*elliptic);
  else if(settings.compareSetting("PRECONDITIONER", "NONE"))
    elliptic->precon = new IdentityPrecon(mesh.Np*mesh.Nelements);

  return *elliptic;
}

elliptic_t::~elliptic_t() {
  maskKernel.free();
  partialAxKernel.free();
  partialGradientKernel.free();
  partialIpdgKernel.free();

  if (traceHalo) traceHalo->Free();
  if (ogsMasked) ogsMasked->Free();
  if (precon) delete precon;
}
