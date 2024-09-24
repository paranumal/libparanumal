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

void cns_t::applyArtificialDiffusion(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T){

// Apply detector which returns an indicator scaled to 0<= s <=1.0
applyDetect(o_Q, o_gradQ, T); 

maxVelocityKernel(mesh.Nelements,
                    BCStateID, 
                     mesh.o_vgeo,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_EToB,
                     o_pCoeff, 
                     o_flowStates, 
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_Q,
                     o_Vmax);


computeViscosityKernel(mesh.Nelements*mesh.Nverts, 
                         avisAlpha, 
                         o_qdetect,
                         o_elmLength,
                         o_Vmax,
                         o_vertexViscosity); 

  
  ogs.GatherScatter(o_vertexViscosity, 1, ogs::Max, ogs::Sym); 
  
  projectViscosityKernel(mesh.Nelements, 
                         o_PM12N,
                         o_vertexViscosity,
                         o_viscosity); 
}


void cns_t::setupArtificialDiffusion(){
  
 // artificial viscosity scale factor:Alpha
  settings.getSetting("ARTDIFF SCALE FACTOR", avisAlpha);  
 
 // Scale with polynomial order
  // avisAlpha = avisAlpha/( (double) std::pow(mesh.N,mesh.dim)); 
  avisAlpha = avisAlpha/( (double) std::pow(mesh.N,1)); 
  
  // Definition of element length might change
  elmLength.malloc(mesh.Nelements); 
  for(dlong e=0; e < mesh.Nelements; e++){
    // printf("%d %.4f %.4f", e, mesh.x[e*mesh.Np + n],mesh.x[e*mesh.Np + n]);
    switch (mesh.elementType) {
      case Mesh::TRIANGLES:
        elmLength[e] = ElementViscosityScaleTri2D(e);
        break;
      case Mesh::QUADRILATERALS:
        elmLength[e] = ElementViscosityScaleQuad2D(e);
        break;
      case Mesh::TETRAHEDRA:
        elmLength[e] = ElementViscosityScaleTet3D(e);
        break;
      case Mesh::HEXAHEDRA:
        elmLength[e] = ElementViscosityScaleHex3D(e);
      break;
      }
    } 
    
  // move to devce
  o_elmLength = platform.malloc<dfloat>(elmLength); 

  // Smooth Viscosity Using Vertex Values
  vertexViscosity.malloc(mesh.Nelements*mesh.Nverts, 0.0); 
  o_vertexViscosity = platform.malloc<dfloat>(vertexViscosity); 
  
  // Allocate Memory for Artificial Viscosity
  viscosity.malloc(mesh.Nelements*mesh.Np, 0.0); 
  o_viscosity = platform.malloc<dfloat>(viscosity); 

  // Allocate Memory for Artificial Viscosity
  Vmax.malloc(mesh.Nelements, 0.0); 
  o_Vmax = platform.malloc<dfloat>(Vmax); 

  memory<dfloat> V, invV1, r1, s1, t1;
  // Compute projection matrix 
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      mesh.NodesTri2D(1, r1, s1);
      mesh.VandermondeTri2D(1, r1, s1, invV1);
      mesh.VandermondeTri2D(1, mesh.r, mesh.s, V);
      break; 
    case Mesh::QUADRILATERALS:
      mesh.NodesQuad2D(1, r1, s1);
      mesh.VandermondeQuad2D(1, r1, s1, invV1);
      mesh.VandermondeQuad2D(1, mesh.r, mesh.s, V);
      break; 
    case Mesh::TETRAHEDRA:
      mesh.NodesTet3D(1, r1, s1, t1);
      mesh.VandermondeTet3D(1, r1, s1, t1, invV1);
      mesh.VandermondeTet3D(1, mesh.r, mesh.s, mesh.t, V);
      break; 
    case Mesh::HEXAHEDRA:
      mesh.NodesHex3D(1,r1, s1, t1); 
      mesh.VandermondeHex3D(1, mesh.r, mesh.s, mesh.t, invV1);
      break; 
  }

  // invert N=1 Vandermonde Matrix
  linAlg_t::matrixInverse(mesh.Nverts, invV1);

  PM12N.malloc(mesh.Nverts*mesh.Np, 0.0); 
  // // Transponse of Projection Operator
  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Nverts; j++){
      dfloat sum = 0; 
      for(int m=0; m<mesh.Nverts; m++){
          sum += V[i*mesh.Nverts + m]*invV1[m*mesh.Nverts + j]; 
      }
      PM12N[j*mesh.Np + i] = sum; 
    }
  } 
  o_PM12N = platform.malloc<dfloat>(PM12N);
  

  // Build gather scatter for vertex-based smmothing
  mesh_t meshC = mesh.SetupNewDegree(1);

  // Define wights for vertex smoothing
  bool verbose = false,  unique  = true; 
  ogs.Setup(mesh.Nelements*mesh.Nverts, meshC.globalIds, mesh.comm, 
                ogs::Signed, ogs::Auto, unique, verbose, platform); 
  
  weight.malloc(mesh.Nelements*mesh.Nverts, 1.0);
  ogs.GatherScatter(weight, 1, ogs::Add, ogs::Sym); 
  // invert weights for gs for vertex nodes only
  for(int i=0; i < mesh.Nelements*mesh.Nverts; i++ ){ 
    weight[i] = 1./weight[i];
  }
  o_weight = platform.malloc<dfloat>(weight); 

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklStabFilePrefix = DCNS "/okl/stab/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  fileName      = oklStabFilePrefix + "artificialDiffusion" + oklFileSuffix;
  kernelName    = "computeViscosity";
  computeViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

  kernelName    = "projectViscosity";
  projectViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

  fileName      = oklStabFilePrefix + "maxVelocity" + suffix + oklFileSuffix;
  kernelName    = "maxVelocity" + suffix;
  maxVelocityKernel  = platform.buildKernel(fileName, kernelName, props);

  if(EulerSolve){
    // kernels from volume file Add isothermal version as well AK. 
    fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
    kernelName = "cnsPartialGradVolumeConservative" + suffix; 
    gradVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
    kernelName = "cnsPartialGradSurfaceConservative" + suffix; // gradient of all conservative fields
    gradSurfaceKernel = platform.buildKernel(fileName, kernelName, props);
}else{
    // kernels from volume file Add isothermal version as well AK. 
    fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradVolumeConservative" + suffix; 
    gradVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradSurfaceConservative" + suffix; // gradient of all conservative fields
    gradSurfaceKernel = platform.buildKernel(fileName, kernelName, props);
 }
 

  if(cubature){    
    if(settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
      if(EulerSolve){
        // kernels from volume file
        fileName   = oklFilePrefix + "cnsCubatureVolumeArtDiff" + suffix + oklFileSuffix;
        kernelName = "cnsEulerCubatureVolumeArtificialDiffsuionLaplace" + suffix;
        cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName, props);
        // kernels from surface file
        fileName   = oklFilePrefix + "cnsCubatureSurfaceArtDiff" + suffix + oklFileSuffix;
        kernelName = "cnsEulerCubatureSurfaceArtificialDiffusionLaplace" + suffix;
        cubatureSurfaceKernel =  platform.buildKernel(fileName, kernelName, props);   
      }else{
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsCubatureVolumeArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsCubatureVolumeArtificialDiffsuionLaplace" + suffix;
      cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

      // kernels from surface file
      fileName   = oklFilePrefix + "cnsCubatureSurfaceArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsCubatureSurfaceArtificialDiffusionLaplace" + suffix;
      cubatureSurfaceKernel =  platform.buildKernel(fileName, kernelName, props);        
      }
    }
   }else{
      if(settings.compareSetting("ARTDIFF TYPE", "LAPLACE")){
      // kernels from volume file
      fileName   = oklFilePrefix + "cnsVolumeArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsVolumeArtificialDiffsuionLaplace" + suffix;
      volumeKernel =  platform.buildKernel(fileName, kernelName, props);


      // kernels from surface file
      fileName   = oklFilePrefix + "cnsSurfaceArtDiff" + suffix + oklFileSuffix;
      kernelName = "cnsSurfaceArtificialDiffusionLaplace" + suffix;
      surfaceKernel =  platform.buildKernel(fileName, kernelName, props);

    }else if(settings.compareSetting("ARTDIFF TYPE", "PHYSICAL")){
      // // kernels from volume file
      // fileName   = oklFilePrefix + "cnsVolumeArtDiff" + suffix + oklFileSuffix;
      // kernelName = "cnsVolumeArtificialDiffsuionPhysical" + suffix;
      // volumeKernel =  platform.buildKernel(fileName, kernelName, props);

      // // kernels from surface file
      // fileName   = oklFilePrefix + "cnsSurfaceArtDiff" + suffix + oklFileSuffix;
      // kernelName = "cnsSurfaceArtificialDiffusionPhysical" + suffix;
      // surfaceKernel =  platform.buildKernel(fileName, kernelName, props);
    }
  }



  // report forces on selected BCs
  fileName   = oklFilePrefix + "cnsForces" + suffix + oklFileSuffix;
  kernelName = "cnsForcesVolume" + suffix;
  forcesVolumeKernel = platform.buildKernel(fileName, kernelName, props);

  kernelName = "cnsForcesSurface" + suffix;
  forcesSurfaceKernel = platform.buildKernel(fileName, kernelName, props);

  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName, props);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    kernelName = "cnsInitialCondition3D";
  }
    initialConditionKernel = platform.buildKernel(fileName, kernelName, props);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "cnsMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, props);

}

 //   // kernels from volume file Add isothermal version as well AK. 
  // fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
  // kernelName = "cnsGradCubatureVolumeConservative" + suffix; 
  // cubatureGradVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

  // // kernels from surface file
  // fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
  // kernelName = "cnsGradCubatureSurfaceConservative" + suffix; // gradient of all conservative fields
  // cubatureGradSurfaceKernel = platform.buildKernel(fileName, kernelName, props);


// Compute h/N for all element types
dfloat cns_t::ElementViscosityScaleTri2D(dlong e) {
  dfloat he = std::numeric_limits<dfloat>::max();
  for(int f=0;f<mesh.Nfaces;++f){
    dlong sid   = mesh.Nsgeo*(mesh.Nfaces*e + f);
    dfloat sJ   = mesh.sgeo[sid + mesh.SJID];
    dfloat invJ = mesh.sgeo[sid + mesh.IJID];
    // sJ = L/2, J = A/2; sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h  = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);
    he = std::min(he, hest);
  }
  return he;
}

dfloat cns_t::ElementViscosityScaleQuad2D(dlong e) {
  dfloat he = std::numeric_limits<dfloat>::max();
  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<mesh.Np;n++){
    J += mesh.vgeo[mesh.Nvgeo*mesh.Np*e + n + mesh.Np*mesh.JWID];
  }

  for(int f=0;f<mesh.Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<mesh.Nfp;i++){
      sJ += mesh.sgeo[mesh.Nsgeo*(mesh.Nfaces*mesh.Nfp*e + mesh.Nfp*f + i) + mesh.WSJID];
    }

    dfloat hest = J/sJ;
    he = std::min(he, hest);
  }
  return he;
}

dfloat cns_t::ElementViscosityScaleTet3D(dlong e) {

  dfloat he = std::numeric_limits<dfloat>::max();
  for(int f=0;f<mesh.Nfaces;++f){
    dlong sid = mesh.Nsgeo*(mesh.Nfaces*e + f);
    dfloat sJ   = mesh.sgeo[sid + mesh.SJID];
    dfloat invJ = mesh.sgeo[sid + mesh.IJID];
    // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);
    he = std::min(he, hest);
  }
  return he;
}

dfloat cns_t::ElementViscosityScaleHex3D(dlong e) {
  dfloat he = std::numeric_limits<dfloat>::max();
  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<mesh.Np;n++){
    J += mesh.vgeo[mesh.Nvgeo*mesh.Np*e + n + mesh.Np*mesh.JWID];
  }
  
  for(int f=0;f<mesh.Nfaces;++f){
      //sum weighted surface Jacobians to integrate over face
      dfloat sJ = 0.0;
      for (int i=0;i<mesh.Nfp;i++){
        sJ += mesh.sgeo[mesh.Nsgeo*(mesh.Nfaces*mesh.Nfp*e + mesh.Nfp*f + i) + mesh.WSJID];
      }

      // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
      // h = 1/(sJ/J)
      dfloat hest = J/sJ;

      he = std::min(he, hest);
  }
  return he;
};






void cns_t::setupLimiter(){
 
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
  props["defines/" "p_Ngrads"]= Ngrads;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  if(cubature){
    // kernels from volume file Add isothermal version as well AK. 
    fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradCubatureVolumeConservative" + suffix; 
    cubatureGradVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradCubatureSurfaceConservative" + suffix; // gradient of all conservative fields
    cubatureGradSurfaceKernel = platform.buildKernel(fileName, kernelName, props);

    // kernels from volume file
    fileName   = oklFilePrefix + "cnsCubatureVolumeLimiter" + suffix + oklFileSuffix;
    kernelName = "cnsCubatureVolumeLimiter" + suffix;
    cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsCubatureSurfaceLimiter" + suffix + oklFileSuffix;
    kernelName = "cnsCubatureSurfaceLimiter" + suffix;
    cubatureSurfaceKernel =  platform.buildKernel(fileName, kernelName, props);
   }else{
    // kernels from volume file Add isothermal version as well AK. 
    fileName   = oklFilePrefix + "cnsGradVolumeConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradVolumeConservative" + suffix; 
    gradVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsGradSurfaceConservative" + suffix + oklFileSuffix;
    kernelName = "cnsGradSurfaceConservative" + suffix; // gradient of all conservative fields
    gradSurfaceKernel = platform.buildKernel(fileName, kernelName, props);

    // kernels from volume file
    fileName   = oklFilePrefix + "cnsVolumeLimiter" + suffix + oklFileSuffix;
    kernelName = "cnsVolumeLimiter" + suffix;
    volumeKernel =  platform.buildKernel(fileName, kernelName, props);

    // kernels from surface file
    fileName   = oklFilePrefix + "cnsSurface" + suffix + oklFileSuffix;
    kernelName = "cnsSurface" + suffix;
    surfaceKernel =  platform.buildKernel(fileName, kernelName, props);
  }

  // Limiter based kernels
  // // projection kernel
  fileName   = oklFilePrefix + "cnsLimiter" + suffix + oklFileSuffix;
  kernelName = "cnsLimiterReconstruct" + suffix;
  limiterReconstructKernel = platform.buildKernel(fileName, kernelName, props);


  kernelName = "cnsLimiterGradient" + suffix;
  limiterGradientKernel = platform.buildKernel(fileName, kernelName, props);


  kernelName = "cnsLimiterVertexBoundary" + suffix;
  limiterVertexBoundaryKernel = platform.buildKernel(fileName, kernelName, props);


  kernelName = "cnsReportArrangeLayout";
  reportArrangeLayoutKernel = platform.buildKernel(fileName, kernelName, props);

  kernelName = "cnsReportAverage";
  reportAverageKernel = platform.buildKernel(fileName, kernelName, props);

  // Report forces on wall bc's
  fileName   = oklFilePrefix + "cnsForces" + suffix + oklFileSuffix;
  kernelName = "cnsForcesVolume" + suffix;
  forcesVolumeKernel = platform.buildKernel(fileName, kernelName, props);

  kernelName = "cnsForcesSurface" + suffix;
  forcesSurfaceKernel = platform.buildKernel(fileName, kernelName, props);

  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName, props);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    kernelName = "cnsInitialCondition3D";
  }
  
  initialConditionKernel = platform.buildKernel(fileName, kernelName, props);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "cnsMaxWaveSpeed" + suffix;
  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, props);

}


void cns_t::setupNoStab(){
  // Needed gradients for velocity only 
  Ngrads = mesh.dim*mesh.dim;

  // setup trace halo exchange 
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  // Update Ngrads for artificial diffusion on device
  props["defines/" "p_Ngrads"]= Ngrads;

  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();

  std::string oklFilePrefix = DCNS "/okl/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  // kernels from volume file Add isothermal version as well AK. 
  fileName   = oklFilePrefix + "cnsGradVolume" + suffix + oklFileSuffix;
  kernelName = "cnsGradVolume" + suffix; 
  gradVolumeKernel =  platform.buildKernel(fileName, kernelName, props);

  // kernels from surface file
  fileName   = oklFilePrefix + "cnsGradSurface" + suffix + oklFileSuffix;
  kernelName = "cnsGradSurface" + suffix; // gradient of all conservative fields
  gradSurfaceKernel = platform.buildKernel(fileName, kernelName, props);

  if(cubature){
    // kernels from volume file
    fileName   = oklFilePrefix + "cnsCubatureVolume" + suffix + oklFileSuffix;
    kernelName = "cnsCubatureVolume" + suffix;
    cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName, props);


    // kernels from surface file
    fileName   = oklFilePrefix + "cnsCubatureSurface" + suffix + oklFileSuffix;
    kernelName = "cnsCubatureSurface" + suffix;
    cubatureSurfaceKernel =  platform.buildKernel(fileName, kernelName, props);
  }else{

     // kernels from volume file
    fileName   = oklFilePrefix + "cnsVolume" + suffix + oklFileSuffix;
    kernelName = "cnsVolume" + suffix;
    volumeKernel =  platform.buildKernel(fileName, kernelName, props);


    // kernels from surface file
    fileName   = oklFilePrefix + "cnsSurface" + suffix + oklFileSuffix;
    kernelName = "cnsSurface" + suffix;
    surfaceKernel =  platform.buildKernel(fileName, kernelName, props);

    
  }


  // vorticity calculation
  fileName   = oklFilePrefix + "cnsVorticity" + suffix + oklFileSuffix;
  kernelName = "cnsVorticity" + suffix;

  vorticityKernel = platform.buildKernel(fileName, kernelName,
                                 props);

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "cnsInitialCondition2D" + oklFileSuffix;
    kernelName = "cnsInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "cnsInitialCondition3D" + oklFileSuffix;
    kernelName = "cnsInitialCondition3D";
  }
    initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                        props);

  fileName   = oklFilePrefix + "cnsMaxWaveSpeed" + suffix + oklFileSuffix;
  kernelName = "cnsMaxWaveSpeed" + suffix;

  maxWaveSpeedKernel = platform.buildKernel(fileName, kernelName, props);

}



// void cns_t::setStabTypes(const stab::detector _Detector, 
//                          const stab::type _Type) {

// // Set Solver Dependent Parameters
//   if (_Solver==Stab::HJS) {
//     solver = Stab::HJS;
//     // Local Hamilton Jacobi Solver with 2 fields i.e. qp and qm;  
//     Ndfields = 2;        // number of fields to be detected
//     Nsfields = 2;        // number of fields to be stabilized
//   } else if (_Solver==Stab::CNS) {
//     solver = Stab::CNS;
//     Ndfields  = 1; 
//     Nsfields  = mesh.dim==2 ? 4:5; // non-isothermal flow only
//   } else if (_Solver==Stab::INS) {
//     LIBP_FORCE_ABORT("Unknown solver type: " << solver);
//   } else {
//     LIBP_FORCE_ABORT("Unknown solver type: " << solver);
//   }

// // Set Detector Type
// if (_Detector==Stab::NODETECT) {
//     detector= Stab::NODETECT;
//  }else if (_Detector==Stab::ALL) {
//     detector = Stab::ALL;
//  }else if (_Detector==Stab::KLOCKNER) {
//     detector = Stab::KLOCKNER;
//  }else if (_Detector==Stab::PERSSON) {
//     detector = Stab::PERSSON;
//  }else if (_Detector==Stab::DUCROS) {
//     detector = Stab::DUCROS;
//  }else {
//     LIBP_FORCE_ABORT("Unknown detector type: " << detector);
//   }


// // Set Stabilization Type
// if (_Type==Stab::FILTER) {
//     type = Stab::FILTER;
//  }else if (_Type==Stab::LIMITER) {
//      type = Stab::LIMITER;
//  }else if (_Type==Stab::ARTDIFF) {
//      type = Stab::ARTDIFF;
//  }else if (_Type==Stab::SUBCELL) {
//      type = Stab::SUBCELL;
//  }else if (_Type==Stab::NOSTAB) {
//      type = Stab::NOSTAB;
//  }else {
//     LIBP_FORCE_ABORT("Unknown solver type: " << type);
//   }

// props["defines/" "s_Ndfields"]= Ndfields;
// props["defines/" "s_Nsfields"]= Nsfields;

// }


// } //namespace libp
