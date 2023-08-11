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

template <typename T>
void printMatrixLocal(linAlgMatrix_t<T> &A, const char *str){

#if 1
  std::cout << "matrix: " <<  std::string(str) << "[" << std::endl;
  for(int r=1;r<=A.rows();++r){
    for(int c=1;c<=A.cols();++c){
      //      std::cout << A(r,c) << " ";
      printf("% 5.4e ", A(r,c));
      
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << "]" << std::endl;
#endif
}


void wave_t::Setup(platform_t& _platform,
                   mesh_t& _mesh,
                   waveSettings_t& _settings){
  
  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  platform.linAlg().InitKernels({"set", "sum", "norm2", "min", "max", "add"});

  
  // FOR THIS - WE ASSUME IPDG  and LAMBDA=1
  Nfields = 1;

  ellipticSettings = _settings.extractEllipticSettings();

  int NBCTypes = 7;
  memory<int> BCType(NBCTypes);
  // bc=3 => outflow => Neumann   => vBCType[3] = 2, etc.
  BCType[0] = 0;
  BCType[1] = 1;
  BCType[2] = 1;
  BCType[3] = 2;
  BCType[4] = 1;
  BCType[5] = 2;
  BCType[6] = 2;

  lambdaSolve = 1;
  elliptic.Setup(platform, mesh, ellipticSettings, lambdaSolve, NBCTypes, BCType);

  // find out if this is a C0 discretization
  disc_c0 = elliptic.settings.compareSetting("DISCRETIZATION", "CONTINUOUS") ? 1 : 0;

  if(disc_c0==1 && mesh.elementType==Mesh::TRIANGLES){
    std::cout << "TRYING TO USE TRIANGLE MESH WITH C0 NOT ALLOWED WITH WAVE" << std::endl;
    exit(-1);
  }


  if(disc_c0==1 && mesh.elementType==Mesh::TETRAHEDRA){
    std::cout << "TRYING TO USE TETRAHEDRA MESH WITH C0 NOT ALLOWED WITH WAVE" << std::endl;
    exit(-1);
  }

  
  // initialize time stepper
  linAlgMatrix_t<dfloat> BTABLE;

  Nstages = 0;
  embedded = 0;

  std::string timeIntegrator;
  settings.getSetting("TIME INTEGRATOR", timeIntegrator);
  libp::TimeStepper::butcherTables(timeIntegrator.c_str(), Nstages, embedded, BTABLE);
  
  // extract alphas, betas
  alpha.reshape(Nstages+1, Nstages);
  beta.reshape(1,Nstages);
  betahat.reshape(1,Nstages);
  esdirkC.reshape(1,Nstages);
  
  for(int n=1;n<=Nstages;++n){
    for(int m=1;m<=Nstages;++m){
      alpha(n,m) = BTABLE(n,m+1);
    }
    beta(1,n) = BTABLE(Nstages+1,n+1);

    if(embedded)
       betahat(1,n) = BTABLE(Nstages+2,n+1);

    esdirkC(n) = BTABLE(n,1);
  }

  for(int m=1;m<=Nstages;++m){
    alpha(Nstages+1,m) = beta(1,m);
  }

  gamma = alpha(2,2);
  invGamma = 1./gamma;
  
  std::cout << "gamma = " << gamma << std::endl;

  alphatilde.reshape(Nstages+1, Nstages);
  gammatilde.reshape(1,Nstages);
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

  for(int i=1;i<=Nstages;++i){
    alphatilde(1+Nstages, i) = betaAlpha(1,i);
  }
  

  
#if 1
  printMatrixLocal(BTABLE, "BUTCHER TABLE");
  printMatrixLocal(alpha, "ALPHA BLOCK");
  printMatrixLocal(beta, "BETA BLOCK");
  printMatrixLocal(gammatilde, "GAMMA TILDE BLOCK");
  printMatrixLocal(alphatilde, "ALPHA TILDE BLOCK");
  printMatrixLocal(betaAlpha,  "BETAxALPHA  BLOCK");
  printMatrixLocal(esdirkC,    "ESDIRK-C    BLOCK");
#endif
  
  o_alphatilde = platform.malloc<dfloat>((Nstages+1)*Nstages, alphatilde.data);
  o_gammatilde = platform.malloc<dfloat>(Nstages, gammatilde.data);
  o_betaAlpha  = platform.malloc<dfloat>(Nstages, betaAlpha.data);
  o_betahatAlpha  = platform.malloc<dfloat>(Nstages, betahatAlpha.data);
  o_esdirkC = platform.malloc<dfloat>(Nstages, esdirkC.data);
  
  o_alpha = platform.malloc<dfloat>((Nstages+1)*Nstages, alpha.data);
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

  kernelInfo["defines/p_Nmax"] = (int) std::max(mesh.Np, mesh.Nfp*mesh.Nfaces);
  
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

  std::string oklFilePrefix = DWAVE "/okl/";
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

#if 0
  kernelName = "waveStepInitializeV2";
  waveStepInitializeKernelV2 =
     platform.buildKernel(fileName, kernelName, kernelInfo);  
  
  kernelName = "waveStageInitializeV2";
  waveStageInitializeKernelV2 =
     platform.buildKernel(fileName, kernelName, kernelInfo);
  
  kernelName = "waveStageFinalizeV2";
  waveStageFinalizeKernelV2 =
     platform.buildKernel(fileName, kernelName, kernelInfo);
#endif  
  
  // element type specific kernels
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

  kernelName = "waveForcing" + suffix;
  waveForcingKernel
     = platform.buildKernel(fileName, kernelName, kernelInfo);


  if(settings.compareSetting("ENABLE FLUX SOURCE","TRUE")){
    // new file
    fileName   = oklFilePrefix + "waveSurfaceSource"  + suffix + oklFileSuffix;
    
    kernelName = "waveSurfaceSource" + suffix;
    std::cout << "kernelName=" << kernelName << std::endl;
    
    waveSurfaceSourceKernel =
       platform.buildKernel(fileName, kernelName, kernelInfo);
  }

  
  //setup linear solver
  Nall = mesh.Np*(mesh.Nelements+mesh.totalHaloPairs);
  
  //create occa buffers
  DL.malloc(Nall);
  PL.malloc(Nall);
  DrhsL.malloc(Nall);
  PhatL.malloc(Nall*Nstages);
  DhatL.malloc(Nall*Nstages);

  o_DL        = platform.malloc<dfloat>(Nall);
  o_PL        = platform.malloc<dfloat>(Nall);
  o_DtildeL   = platform.malloc<dfloat>(Nall);
  o_DrhsL     = platform.malloc<dfloat>(Nall);
  o_DhatL     = platform.malloc<dfloat>(Nall*Nstages);
  o_PhatL     = platform.malloc<dfloat>(Nall*Nstages);
  o_scratch1L = platform.malloc<dfloat>(Nall);
  o_scratch2L = platform.malloc<dfloat>(Nall);

  o_FL        = platform.malloc<dfloat>(Nall);
  o_filtPL    = platform.malloc<dfloat>(Nall);
  
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
    invMM.malloc(mesh.Np*mesh.Np);
    MM.malloc(mesh.Np*mesh.Np);
    int cnt = 0;
    for(int j=0;j<mesh.Nq;++j){
      for(int i=0;i<mesh.Nq;++i){
	MM[cnt] = mesh.gllw[i]*mesh.gllw[j];
	invMM[cnt] = 1./MM[cnt];
	++cnt;
      }
    }
      
  } else if(mesh.elementType==Mesh::TETRAHEDRA) {
    mesh.VandermondeTet3D(mesh.N, mesh.r, mesh.s, mesh.t, V);
    mesh.invMassMatrixTet3D(mesh.Np, V, invMM);
    mesh.MassMatrixTet3D(mesh.Np, V, MM);
  } else {
    invMM.malloc(mesh.Np*mesh.Np);
    MM.malloc(mesh.Np*mesh.Np);

    int cnt = 0;
    for(int k=0;k<mesh.Nq;++k){
      for(int j=0;j<mesh.Nq;++j){
	for(int i=0;i<mesh.Nq;++i){
	  MM[cnt] = mesh.gllw[i]*mesh.gllw[j]*mesh.gllw[k];
	  invMM[cnt] = 1./MM[cnt];
	  printf("MM[%d]=%g\n", cnt, MM[cnt]);
	  ++cnt;
	}
      }
    }
    
    
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


  dfloat hmin = 1e9, hmax = -1e9;
  for(dlong e=0;e<mesh.Nelements;++e){
    dfloat h = mesh.ElementCharacteristicLength(e);
    hmin = std::min(hmin, h);
    hmax = std::max(hmax, h);
  }

  std::cout << "h in range [ " << hmin << ", " << hmax << "] " << std::endl;

  if(settings.compareSetting("ENABLE FLUX SOURCE","TRUE")){
  /* build source via fluxes */
  dfloat rdim = 1;
  xsource = 0;
  ysource = 0;
  zsource = 0;
  fsource = 1;
  settings.getSetting("FLUX SOURCE X", xsource);
  settings.getSetting("FLUX SOURCE Y", ysource);
  settings.getSetting("FLUX SOURCE Z", zsource);
  settings.getSetting("FLUX SOURCE PATCH RADIUS", rdim);

  settings.getSetting("FLUX SOURCE FREQUENCY", fsource);

  patchLabels.malloc(mesh.Nelements);
  for(dlong e1=0;e1<mesh.Nelements;++e1){
    patchLabels[e1] = 0;
  }

  for(dlong e=0;e<mesh.Nelements;++e){

    // test if source is in this element
    if(mesh.PointInclusionTest(e, xsource, ysource, 0)){
      std::cout << "Found element " << e <<
         " includes source ("  << xsource << "," << ysource << ")" << std::endl;
      patchLabels[e] = 1;
#if 1
      for(int f1=0;f1<mesh.Nfaces;++f1){
        dlong e1 = mesh.EToE[e*mesh.Nfaces+f1];
        patchLabels[e1] = 1;
        for(int f2=0;f2<mesh.Nfaces;++f2){
          dlong e2 = mesh.EToE[e1*mesh.Nfaces+f2];
          patchLabels[e2] = 1;
          for(int f3=0;f3<mesh.Nfaces;++f3){
            dlong e3 = mesh.EToE[e2*mesh.Nfaces+f3];
            patchLabels[e3] = 1;
          }
        }
      }
#endif
    }
  }

  int change = 1;
  while(change>0){
    change = 0;
    for(dlong e=0;e<mesh.Nelements;++e){
      if(patchLabels[e]==0){
        int countNeighbors = 0;
        for(int f=0;f<mesh.Nfaces;++f){
          dlong eP = mesh.EToE[e*mesh.Nfaces+f];
          if(patchLabels[eP]==1)
             ++countNeighbors;
        }
        if(countNeighbors==2){
           patchLabels[e] = 1;
           ++change;
        }
      }
    }
  }
  memory<int> EToPatch(mesh.Nelements*mesh.Nfaces);
  for(dlong e1=0;e1<mesh.Nelements;++e1){
    for(int f1=0;f1<mesh.Nfaces;++f1){
      EToPatch[e1*mesh.Nfaces+f1] = 0;
    }
  }
  
  for(dlong e1=0;e1<mesh.Nelements;++e1){
    for(int f1=0;f1<mesh.Nfaces;++f1){
      dlong e2 = mesh.EToE[e1*mesh.Nfaces+f1];
      
      // only handle internal patches first
      if(patchLabels[e1]==1){
        if(patchLabels[e2]==0){
          EToPatch[e1*mesh.Nfaces+f1] = 1;
        }
      }
      if(patchLabels[e1]==0){
        if(patchLabels[e2]==1){
          EToPatch[e1*mesh.Nfaces+f1] = -1;
        }
      }
    }
  }

  for(dlong e1=0;e1<mesh.Nelements;++e1){
    for(int f1=0;f1<mesh.Nfaces;++f1){
      if(EToPatch[e1*mesh.Nfaces+f1]){
        dlong e2 = mesh.EToE[e1*mesh.Nfaces+f1];
        dlong f2 = mesh.EToF[e1*mesh.Nfaces+f1];
        printf("(%d,%d) patch %d =>  neighbor: ( %d,%d) patch %d\n",
               e1, f1, EToPatch[e1*mesh.Nfaces+f1],
               e2, f2, EToPatch[e2*mesh.Nfaces+f2]);
      }
    }
  }

  o_EToPatch = platform.malloc<int>(mesh.Nelements*mesh.Nfaces, EToPatch);
  }
}
