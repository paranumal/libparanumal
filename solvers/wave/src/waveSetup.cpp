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

#include "wave.hpp"

void wave_t::Setup(platform_t& _platform,
                   mesh_t& _mesh,
                   waveSettings_t& _settings){
  
  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  ellipticSettings = _settings.extractEllipticSettings();
  
  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"set", "sum", "norm2", "min", "max"});
  
  // FOR THIS - WE ASSUME IPDG  and LAMBDA=1
  Nfields = 1;
  
  // initialize time stepper
  linAlgMatrix_t<dfloat> BTABLE;

  Nstages = 0;
  embedded = 0;
  
  libp::TimeStepper::butcherTables("ESDIRK6(5)9L[2]SA", Nstages, embedded, BTABLE);
  
  // extract alphas, betas
  alpha.reshape(Nstages, Nstages);
  beta.reshape(1,Nstages);
  betahat.reshape(1,Nstages);
  
  for(int n=1;n<=Nstages;++n){
    for(int m=1;m<=Nstages;++m){
      alpha(n,m) = BTABLE(n,m+1);
    }
    beta(1,n) = BTABLE(Nstages+1,n+1);

    if(embedded)
       betahat(1,n) = BTABLE(Nstages+2,n+1);

  }

  gamma = alpha(2,2);
  dfloat invGamma = 1./gamma;
  
  std::cout << "gamma = " << gamma << std::endl;

  linAlgMatrix_t<dfloat> alphatilde(Nstages, Nstages);
  linAlgMatrix_t<dfloat> gammatilde(1,Nstages);
  linAlgMatrix_t<dfloat> betaAlpha(1,Nstages);
  linAlgMatrix_t<dfloat> betahatAlpha(1,Nstages);
  alphatilde = 0.;
  gammatilde = 0.;
  
  for(int i=0;i<=Nstages-1;++i){
    for(int j=1;j<=i;++j){
      alphatilde(i+1,j) += alpha(i+1,j)/gamma;
    }
    for(int j=1;j<=i;++j){
      for(int k=1;k<=j;++k){
        alphatilde(i+1,k) += alpha(i+1,j)*alpha(j,k)/(gamma*gamma);
      }
    }
    gammatilde(1,i+1) = 1./gamma;
    for(int j=1;j<=i;++j){
      gammatilde(1,i+1) += alpha(i+1,j)/(gamma*gamma);
    }
  }

  betaAlpha = 0.;
  betahatAlpha = 0.;
  for(int i=1;i<=Nstages;++i){
    for(int j=1;j<=i;++j){
      betaAlpha(1,j) += beta(1,i)*alpha(i,j);
      betahatAlpha(1,j) += betahat(1,i)*alpha(i,j);
    }
  }

#if 0
  printMatrixLocal(BTABLE, "BUTCHER TABLE");
  printMatrixLocal(alpha, "ALPHA BLOCK");
  printMatrixLocal(beta, "BETA BLOCK");
  printMatrixLocal(gammatilde, "GAMMA TILDE BLOCK");
  printMatrixLocal(alphatilde, "ALPHA TILDE BLOCK");
  printMatrixLocal(betaAlpha,  "BETAxALPHA  BLOCK");
#endif
  
  o_alphatilde = platform.malloc<dfloat>(Nstages*Nstages, alphatilde.data);
  o_gammatilde = platform.malloc<dfloat>(Nstages, gammatilde.data);
  o_betaAlpha  = platform.malloc<dfloat>(Nstages, betaAlpha.data);
  o_betahatAlpha  = platform.malloc<dfloat>(Nstages, betahatAlpha.data);
  
  o_alpha = platform.malloc<dfloat>(Nstages*Nstages, alpha.data);
  o_beta  = platform.malloc<dfloat>(Nstages, beta.data);
  o_betahat  = platform.malloc<dfloat>(Nstages, betahat.data);
  
  properties_t kernelInfo = mesh.props; //copy base occa properties

  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  //add standard boundary functions
  std::string boundaryHeaderFileName;
  if (mesh.dim==2)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary2D.h");
  else if (mesh.dim==3)
    boundaryHeaderFileName = std::string(DELLIPTIC "/data/ellipticBoundary3D.h");
  kernelInfo["includes"] += boundaryHeaderFileName;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES) {
    suffix = "Tri2D";
  } else if(mesh.elementType==Mesh::QUADRILATERALS) {
    suffix = "Quad2D";
  } else if(mesh.elementType==Mesh::TETRAHEDRA) {
    suffix = "Tet3D";
  } else { //mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";
  }

  std::string oklFilePrefix = DELLIPTIC "/okl/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  // WE WILL ASSUME ZERO DIRICHLET BOUNDARIES
  
  fileName   = oklFilePrefix + "waveKernels"  + oklFileSuffix;
  
  kernelName = "waveStageUpdate";
  waveStageUpdateKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "waveCombine";
  waveCombineKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "waveErrorEstimate";
  waveErrorEstimateKernel
     = platform.buildKernel(fileName, kernelName, kernelInfo);
  
  fileName   = oklFilePrefix + "waveKernels" + suffix + oklFileSuffix;

  kernelName = "waveStepInitialize" + suffix;
  waveStepInitializeKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "waveStageFinalize" + suffix;
  waveStageFinalizeKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);
  
  kernelName = "waveStepFinalize" + suffix;
  waveStepFinalizeKernel =
     platform.buildKernel(fileName, kernelName, kernelInfo);
  
  kernelName = "waveInitialConditions" + suffix;
  waveInitialConditionsKernel
     = platform.buildKernel(fileName, kernelName, kernelInfo);

  //create occa buffers
  DL.malloc(Nall);
  PL.malloc(Nall);
  DrhsL.malloc(Nall);
  PhatL.malloc(Nall);
  DhatL.malloc(Nall*Nstages);

  o_DL        = platform.malloc<dfloat>(Nall);
  o_PL        = platform.malloc<dfloat>(Nall);
  o_DtildeL   = platform.malloc<dfloat>(Nall);
  o_DrhsL     = platform.malloc<dfloat>(Nall);
  o_DhatL     = platform.malloc<dfloat>(Nall*Nstages);
  o_PhatL     = platform.malloc<dfloat>(Nall);
  o_scratch1L = platform.malloc<dfloat>(Nall);
  o_scratch2L = platform.malloc<dfloat>(Nall);

  if (disc_c0){
    NglobalDofs = elliptic.ogsMasked.NgatherGlobal*Nfields;
  } else {
    NglobalDofs = mesh.NelementsGlobal*mesh.Np*Nfields;
  }
  
  o_Dtilde   = platform.malloc<dfloat>(NglobalDofs);
  o_Drhs     = platform.malloc<dfloat>(NglobalDofs);
  o_scratch1 = platform.malloc<dfloat>(NglobalDofs);
  o_scratch2 = platform.malloc<dfloat>(NglobalDofs);
  
  memory<dfloat> invMM, V, MM;
  
  if(mesh.elementType==Mesh::TRIANGLES) {
    mesh.VandermondeTri2D(mesh.N, mesh.r, mesh.s, V);
    mesh.invMassMatrixTri2D(mesh.Np, V, invMM);
    mesh.MassMatrixTri2D(mesh.Np, V, MM);
  } else if(mesh.elementType==Mesh::QUADRILATERALS) {
    invMM.malloc(1);
    MM.malloc(1);
  } else if(mesh.elementType==Mesh::TETRAHEDRA) {
    mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
    mesh.invMassMatrixTet3D(mesh.Np, V, invMM);
    mesh.MassMatrixTet3D(mesh.Np, V, MM);
  } else { 
    printf("HEXES NOT IMPLEMENTED YET\n");
    exit(-1);
  }
  
  o_invMM = platform.malloc<dfloat>(mesh.Np*mesh.Np, invMM);
  o_MM    = platform.malloc<dfloat>(mesh.Np*mesh.Np, MM);

  // triangle specific
  if(mesh.elementType==Mesh::TRIANGLES ||
     mesh.elementType==Mesh::TETRAHEDRA) {

    WJ.malloc(mesh.Nelements);
    invWJ.malloc(mesh.Nelements);

    for(int e=0;e<mesh.Nelements;++e){
      invWJ[e] = 1./mesh.wJ[e];
      WJ[e] = mesh.wJ[e];
    }

    o_invWJ = platform.malloc<dfloat>(mesh.Nelements, invWJ);
    o_WJ    = platform.malloc<dfloat>(mesh.Nelements, WJ);  

  }else{

    // use ungathered weights in WJ on device
    WJ.malloc(mesh.Np*mesh.Nelements);
    for(int n=0;n<mesh.Np*mesh.Nelements;++n){
      WJ[n] = mesh.wJ[n];
    }
    o_WJ    = platform.malloc<dfloat>(mesh.Nelements*mesh.Np, WJ);
    
    // gather weights for C0
    if(disc_c0){
      elliptic.ogsMasked.GatherScatter(WJ, 1, ogs::Add, ogs::Sym);
    }
    
    // use globalized mass for C0
    invWJ.malloc(mesh.Np*mesh.Nelements);
    for(int n=0;n<mesh.Np*mesh.Nelements;++n){
      invWJ[n] = (WJ[n]) ? 1./WJ[n]:0;
    }
    
    o_invWJ = platform.malloc<dfloat>(mesh.Nelements*mesh.Np, invWJ);
  }
}
