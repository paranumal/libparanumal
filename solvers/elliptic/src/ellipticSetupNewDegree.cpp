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

  // elliptic_t* elliptic = new elliptic_t(meshC, linAlg, lambda);
  elliptic_t* elliptic = new elliptic_t(meshC, linAlg);

  //shallow copy
  elliptic->lambda = lambda;

  elliptic->disc_ipdg = disc_ipdg;
  elliptic->disc_c0 = disc_c0;
  elliptic->var_coef = var_coef; 

  elliptic->grad = grad;
  elliptic->o_grad = o_grad;

  elliptic->BCType = BCType;

  elliptic->maskKernel = maskKernel;


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

  // Add coefficients here ......
  if(elliptic->var_coef){
  
 #if 1
 this->SetupNewCoefficient(*elliptic); 
 #else
  // this needs to be changed later
  const dlong Nall   = meshC.Np *(meshC.Nelements+meshC.totalHaloPairs);    

  // currently scalar coefficients are supported
  elliptic->coeff     = (dfloat *) calloc(2*Nall, sizeof(dfloat)); 
  elliptic->o_coeff   = meshC.device.malloc(2*Nall*sizeof(dfloat), elliptic->coeff);

  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  sprintf(fileName, DELLIPTIC "/okl/ellipticCoefficient%s.okl", suffix);
  sprintf(kernelName,"ellipticCoefficient%s", suffix);

    elliptic->coefficientKernel  = buildKernel(meshC.device,fileName, kernelName,
                                              kernelInfo, meshC.comm); 

    elliptic->coefficientKernel(meshC.Nelements,
                                meshC.o_x,
                                meshC.o_y,
                                meshC.o_z,
                                Nall,  
                                elliptic->o_coeff); 
    
    // copy to host for setup

   elliptic->o_coeff.copyTo(elliptic->coeff);
#endif
  }else{ // setting contant coefficient 
  // copy base elliptic coefficient
  elliptic->coeff     = coeff; 
  elliptic->o_coeff   = o_coeff;
  }

  // Ax kernel
  if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
    if(meshC.elementType==HEXAHEDRA){
       if(settings.compareSetting("ELEMENT MAP", "TRILINEAR")){
        if(elliptic->var_coef)
        sprintf(kernelName, "ellipticPartialAxTrilinearVar%s", suffix);
        else
        sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
      }else{
        if(elliptic->var_coef)
        sprintf(kernelName, "ellipticPartialAxVar%s", suffix);
        else
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
      }
    } else{
      if(elliptic->var_coef)      
          sprintf(kernelName, "ellipticPartialAxVar%s", suffix);
        else
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
  
  return *elliptic;
}
