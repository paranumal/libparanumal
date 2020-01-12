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

elliptic_t& elliptic_t::SetupNewDegree(mesh_t& meshC){

  //if asking for the same degree, return the original solver
  if (meshC.N == mesh.N) return *this;

  elliptic_t* elliptic = new elliptic_t(meshC, linAlg, settings, lambda);

  //shallow copy
  elliptic->Nfields = Nfields;
  elliptic->lambda = lambda;

  elliptic->disc_ipdg = disc_ipdg;
  elliptic->disc_c0 = disc_c0;

  elliptic->grad = grad;
  elliptic->o_grad = o_grad;

  elliptic->BCType = BCType;

  elliptic->maskKernel = maskKernel;

  /*setup trace halo exchange */
  elliptic->traceHalo = meshC.HaloTraceSetup(elliptic->Nfields);

  //setup boundary flags and make mask and masked ogs
  elliptic->BoundarySetup();

  //tau (penalty term in IPDG)
  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    if (meshC.elementType==TRIANGLES ||
        meshC.elementType==QUADRILATERALS){
      elliptic->tau = 2.0*(meshC.N+1)*(meshC.N+2)/2.0;
      if(meshC.dim==3)
        elliptic->tau *= 1.5;
    } else
      elliptic->tau = 2.0*(meshC.N+1)*(meshC.N+3);

    //buffer for gradient (Reuse the original buffer)
    // dlong Ntotal = meshC.Np*(meshC.Nelements+meshC.totalHaloPairs);
    // elliptic->grad = (dfloat*) calloc(Ntotal*4, sizeof(dfloat));
    // elliptic->o_grad  = meshC.device.malloc(Ntotal*4*sizeof(dfloat), elliptic->grad);
  }

  // OCCA build stuff
  occa::properties kernelInfo = elliptic->props; //copy base occa properties

  // set kernel name suffix
  char *suffix;
  if(meshC.elementType==TRIANGLES){
    if(meshC.dim==2)
      suffix = strdup("Tri2D");
    else
      suffix = strdup("Tri3D");
  } else if(meshC.elementType==QUADRILATERALS){
    if(meshC.dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  } else if(meshC.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  else if(meshC.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  //add standard boundary functions
  char *boundaryHeaderFileName;
  if (meshC.dim==2)
    boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (meshC.dim==3)
    boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  int NblockV = mymax(1,512/meshC.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  // Ax kernel
  if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
    if(meshC.elementType==HEXAHEDRA){
      if(mesh.settings.compareSetting("ELEMENT MAP", "TRILINEAR"))
        sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
      else
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
    } else{
      sprintf(kernelName, "ellipticPartialAx%s", suffix);
    }

    elliptic->partialAxKernel = buildKernel(meshC.device, fileName, kernelName,
                                            kernelInfo, meshC.comm);

  } else if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    int Nmax = mymax(meshC.Np, meshC.Nfaces*meshC.Nfp);
    kernelInfo["defines/" "p_Nmax"]= Nmax;

    sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
    sprintf(kernelName, "ellipticPartialGradient%s", suffix);
    elliptic->partialGradientKernel = buildKernel(meshC.device, fileName, kernelName,
                                                  kernelInfo, meshC.comm);

    sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
    sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
    elliptic->partialIpdgKernel = buildKernel(meshC.device, fileName, kernelName,
                                              kernelInfo, meshC.comm);
  }

  elliptic->precon = NULL;

  return *elliptic;
}
