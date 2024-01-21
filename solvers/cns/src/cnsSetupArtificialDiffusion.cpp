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

 // kernels from volume file Add isothermal version as well AK. 
  fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
  kernelName = "cnsGradVolumeConservative" + suffix; 
  gradVolumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

  // kernels from surface file
  fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
  kernelName = "cnsGradSurfaceConservative" + suffix; // gradient of all conservative fields
  gradSurfaceKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  if(cubature){
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
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsVolumeArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsVolumeArtificialDiffsuionPhysical" + suffix;
      volumeKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);

      // kernels from surface file
      fileName   = oklFilePrefix + "cnsSurfaceArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsSurfaceArtificialDiffusionPhysical" + suffix;
      surfaceKernel =  platform.buildKernel(fileName, kernelName, kernelInfo);
    }
  }



  // report forces on wall bc's
  fileName   = oklFilePrefix + "cnsComputeForces" + suffix + oklFileSuffix;
  kernelName = "cnsComputeForces" + suffix;
  computeForcesKernel = platform.buildKernel(fileName, kernelName, kernelInfo);

  kernelName = "cnsComputeMoments" + suffix;
  computeMomentsKernel = platform.buildKernel(fileName, kernelName, kernelInfo);



  
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