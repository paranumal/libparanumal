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

#include "stab.hpp"
#include "mesh.hpp"

namespace libp {
// using namespace libp;

void stab_t::stabSetupArtdiff(){
  // Compute viscosity scaling i.e. alpha * mesh_length/ mesh.N^2
  // read viscosity scaling factor i.e. viscosity_max = factor*h/N 
  settings.getSetting("ARTDIFF SCALE FACTOR", scaleFactor);  

  viscosityScale.malloc(mesh.Nelements); 
  for(dlong e=0; e < mesh.Nelements; e++){
    switch (mesh.elementType) {
      case Mesh::TRIANGLES:
        viscosityScale[e] = ElementViscosityScaleTri2D(e);
        break;
      case Mesh::QUADRILATERALS:
        viscosityScale[e] = ElementViscosityScaleQuad2D(e);
        break;
      case Mesh::TETRAHEDRA:
        viscosityScale[e] = ElementViscosityScaleTet3D(e);
        break;
      case Mesh::HEXAHEDRA:
        viscosityScale[e] = ElementViscosityScaleHex3D(e);
        break;
    }
  } 
  o_viscosityScale = platform.malloc<dfloat>(viscosityScale); 

  // Smooth Viscosity Using Vertex Values
  vertexViscosity.malloc(mesh.Nelements*mesh.Nverts*Ndfields, 0.0); 
  o_vertexViscosity = platform.malloc<dfloat>(vertexViscosity); 

  // Allocate Memory for Artificial Viscosity
  viscosity.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*Ndfields, 0.0); 
  o_viscosity = platform.malloc<dfloat>(viscosity); 

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

  projectViscosity.malloc(mesh.Nverts*mesh.Np, 0.0); 
  // // Transponse of Projection Operator
  for(int i=0; i<mesh.Np; i++){
    for(int j=0; j<mesh.Nverts; j++){
      dfloat sum = 0; 
      for(int m=0; m<mesh.Nverts; m++){
          sum += V[i*mesh.Nverts + m]*invV1[m*mesh.Nverts + j]; 
      }
      projectViscosity[j*mesh.Np + i] = sum; 
    }
  } 
  
  o_projectViscosity = platform.malloc<dfloat>(projectViscosity);
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

  
  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int sNblockV = std::max(1, blockMax/mesh.Np);
  props["defines/" "s_NblockV"]= sNblockV;
  props["defines/" "s_Nverts"] = int(mesh.Nverts);

  std::string oklFilePrefix = STAB_DIR "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;
  fileName      = oklFilePrefix + "artificialDiffusion" + oklFileSuffix;
  kernelName    = "computeViscosity";
  computeViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

  kernelName    = "projectViscosity";
  projectViscosityKernel  = platform.buildKernel(fileName, kernelName, props);

  kernelName    = "maskElements";
  maskElementsKernel  = platform.buildKernel(fileName, kernelName, props);
}


void stab_t::stabApplyArtdiff(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // detect elements
  Detect(o_Q, o_RHS, T); 

  computeViscosityKernel(mesh.Nelements*mesh.Nverts, 
                         scaleFactor, 
                         o_viscosityActivation,
                         o_viscosityScale,
                         o_vertexViscosity); 

  // // 
  ogs.GatherScatter(o_vertexViscosity, 1, ogs::Add, ogs::Sym); 
  platform.linAlg().amx(mesh.Nelements*mesh.Nverts, 1.0, o_weight, o_vertexViscosity); 

  projectViscosityKernel(mesh.Nelements, 
                         o_projectViscosity,
                         o_vertexViscosity,
                         o_viscosity); 


}


// Compute h/N for all element types
dfloat stab_t::ElementViscosityScaleTri2D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();
  for(int f=0;f<mesh.Nfaces;++f){
    dlong sid = mesh.Nsgeo*(mesh.Nfaces*e + f);
    dfloat sJ   = mesh.sgeo[sid + mesh.SJID];
    dfloat invJ = mesh.sgeo[sid + mesh.IJID];

    // sJ = L/2, J = A/2; sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h  = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);
    h = std::min(h, hest);
  }
  return h/(mesh.N*mesh.N);
}

dfloat stab_t::ElementViscosityScaleQuad2D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();

  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<mesh.Np;n++)
    J += mesh.vgeo[mesh.Nvgeo*mesh.Np*e + n + mesh.Np*mesh.JWID];

  for(int f=0;f<mesh.Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<mesh.Nfp;i++)
      sJ += mesh.sgeo[mesh.Nsgeo*(mesh.Nfaces*mesh.Nfp*e + mesh.Nfp*f + i) + mesh.WSJID];

    // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
    // h = 1/(sJ/J)
    dfloat hest = J/sJ;

    h = std::min(h, hest);
  }
  return h/(mesh.N*mesh.N);
}

dfloat stab_t::ElementViscosityScaleTet3D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();
  for(int f=0;f<mesh.Nfaces;++f){
    dlong sid = mesh.Nsgeo*(mesh.Nfaces*e + f);
    dfloat sJ   = mesh.sgeo[sid + mesh.SJID];
    dfloat invJ = mesh.sgeo[sid + mesh.IJID];

    // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);

    h = std::min(h, hest);
  }
  return h/(mesh.N*mesh.N*mesh.N);
}

dfloat stab_t::ElementViscosityScaleHex3D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();

  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<mesh.Np;n++)
    J += mesh.vgeo[mesh.Nvgeo*mesh.Np*e + n + mesh.Np*mesh.JWID];

  for(int f=0;f<mesh.Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<mesh.Nfp;i++)
      sJ += mesh.sgeo[mesh.Nsgeo*(mesh.Nfaces*mesh.Nfp*e + mesh.Nfp*f + i) + mesh.WSJID];

    // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
    // h = 1/(sJ/J)
    dfloat hest = J/sJ;

    h = std::min(h, hest);
  }
  return h/(mesh.N*mesh.N*mesh.N);
};


} //namespace libp
