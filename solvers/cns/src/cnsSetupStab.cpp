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

#include "cns.hpp"

void cns_t::setupArtificialDiffusion(properties_t & kernelInfo){
 
  // needed gradients of all fields 
  Ngrads = Nfields*mesh.dim;

  // setup trace halo exchange 
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  // Update Ngrads for artificial diffusion on device
  kernelInfo["defines/" "p_Ngrads"]= Ngrads;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  if(cubature){
    // kernels from volume file Add isothermal version as well AK. 
      fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
      kernelName = "cnsGradCubatureVolumeConservative" + suffix; 
      cubatureGradVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

      // kernels from surface file
      fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
      kernelName = "cnsGradCubatureSurfaceConservative" + suffix; // gradient of all conservative fields
      cubatureGradSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

      if(settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
        // kernels from volume file
        fileName   = oklFilePrefix + "cnsCubatureVolumeArtDiff" + suffix + oklFileSuffix;
        kernelName = "cnsCubatureVolumeArtificialDiffsuionLaplace" + suffix;
        cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

        // kernels from surface file
        fileName   = oklFilePrefix + "cnsCubatureSurfaceArtDiff" + suffix + oklFileSuffix;
        kernelName = "cnsCubatureSurfaceArtificialDiffusionLaplace" + suffix;
        cubatureSurfaceKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);
      }
   }else{
      // kernels from volume file Add isothermal version as well AK. 
      fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
      kernelName = "cnsGradVolumeConservative" + suffix; 
      gradVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

      // kernels from surface file
      fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
      kernelName = "cnsGradSurfaceConservative" + suffix; // gradient of all conservative fields
      gradSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

      if(settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsVolumeArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsVolumeArtificialDiffsuionLaplace" + suffix;
      volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);


      // kernels from surface file
      fileName   = oklFilePrefix + "cnsSurfaceArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsSurfaceArtificialDiffusionLaplace" + suffix;
      surfaceKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

    }else if(settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
      // // kernels from volume file
      // fileName   = oklFilePrefix + "cnsVolumeArtDiff" + suffix + oklFileSuffix;
      // kernelName = "cnsVolumeArtificialDiffsuionPhysical" + suffix;
      // volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

      // // kernels from surface file
      // fileName   = oklFilePrefix + "cnsSurfaceArtDiff" + suffix + oklFileSuffix;
      // kernelName = "cnsSurfaceArtificialDiffusionPhysical" + suffix;
      // surfaceKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);
    }
  }



  // report forces on wall bc's
  fileName   = oklFilePrefix + "cnsForces" + suffix + oklFileSuffix;
  kernelName = "cnsForcesVolume" + suffix;
  forcesVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "cnsForcesSurface" + suffix;
  forcesSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  
  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    kernelName = "cnsInitialCondition3D";
  }
    initialConditionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "cnsMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

}



void cns_t::setupLimiter(properties_t & kernelInfo){
 
  // needed gradients of all fields 
  Ngrads = Nfields*mesh.dim;

  // Define wights for vertex smoothing
  bool verbose = false,  unique  = true; 
  ogs.Setup(mesh.Nelements*mesh.Np, mesh.globalIds, mesh.comm, 
                ogs::Signed, ogs::Auto, unique, verbose, platform); 

  
  weight.malloc(mesh.Nelements*mesh.Np, 1.0);
  ogs.GatherScatter(weight, 1, ogs::Add, ogs::Sym); 

  // invert weights for gs for vertex nodes only
  for(int i=0; i < mesh.Nelements*mesh.Np; i++ ){ 
    weight[i] = 1./weight[i];
  }
  // Move to device
  o_weight = platform.malloc<dfloat>(weight); 

  // setup trace halo exchange 
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);
  //
  // Update Ngrads for artificial diffusion on device
  kernelInfo["defines/" "p_Ngrads"]= Ngrads;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  if(cubature){
    // kernels from volume file Add isothermal version as well AK. 
    fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradCubatureVolumeConservative" + suffix; 
    cubatureGradVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradCubatureSurfaceConservative" + suffix; // gradient of all conservative fields
    cubatureGradSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    // kernels from volume file
    fileName   = oklFilePrefix + "cnsCubatureVolumeLimiter" + suffix + oklFileSuffix;
    kernelName = "cnsCubatureVolumeLimiter" + suffix;
    cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsCubatureSurfaceLimiter" + suffix + oklFileSuffix;
    kernelName = "cnsCubatureSurfaceLimiter" + suffix;
    cubatureSurfaceKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);
   }else{
    // kernels from volume file Add isothermal version as well AK. 
    fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradVolumeConservative" + suffix; 
    gradVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradSurfaceConservative" + suffix; // gradient of all conservative fields
    gradSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

    // kernels from volume file
    fileName   = oklFilePrefix + "cnsVolumeLimiter" + suffix + oklFileSuffix;
    kernelName = "cnsVolumeLimiter" + suffix;
    volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsSurface" + suffix + oklFileSuffix;
    kernelName = "cnsSurface" + suffix;
    surfaceKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);
  }

  // Limiter based kernels
  // // projection kernel
  fileName   = oklFilePrefix + "cnsLimiter" + suffix + oklFileSuffix;
  kernelName = "cnsLimiterReconstruct" + suffix;
  limiterReconstructKernel = platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "cnsLimiterGradient" + suffix;
  limiterGradientKernel = platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "cnsLimiterVertexBoundary" + suffix;
  limiterVertexBoundaryKernel = platform.buildKernel(fileName, kernelName, kernelInfo);


  kernelName = "cnsReportArrangeLayout";
  reportArrangeLayoutKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "cnsReportAverage";
  reportAverageKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  // Report forces on wall bc's
  fileName   = oklFilePrefix + "cnsForces" + suffix + oklFileSuffix;
  kernelName = "cnsForcesVolume" + suffix;
  forcesVolumeKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "cnsForcesSurface" + suffix;
  forcesSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    kernelName = "cnsInitialCondition3D";
  }
  
  initialConditionKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "cnsMaxWaveSpeed" + suffix;
  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

}