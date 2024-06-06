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

void stab_t::stabSetupLimiter(){
  // Setup Operators
  limiterSetupOperators(); 
  // Setup Geometry
   switch (mesh.elementType) {
      case Mesh::TRIANGLES:
        limiterGeometricFactorsTri2D();
        break;
      case Mesh::QUADRILATERALS:
        // limiterGeometricFactorsQuad2D(); 
        LIBP_FORCE_ABORT("Limiter-Quad is not implemented yet");
        break;
      case Mesh::TETRAHEDRA:
        // limiterGeometricFactorsTet3D(); 
        LIBP_FORCE_ABORT("Limiter-Tet is not implemented yet");
        break;
      case Mesh::HEXAHEDRA:
        // limiterGeometricFactorsHex3D();
        LIBP_FORCE_ABORT("Limiter-Hex is not implemented yet");
        break;
    }

    vertexNodes.malloc(mesh.Nverts); 
    vertexNodes.copyFrom(mesh.vertexNodes); 
    o_vertexNodes = platform.malloc<int>(vertexNodes);

    qc.malloc((mesh.Nelements + mesh.totalHaloPairs)*Nsfields); 
    o_qc = platform.malloc<dfloat>(qc); 

    qv.malloc(mesh.Nelements*Nsfields*mesh.Nverts); 
    o_qv = platform.malloc<dfloat>(qv); 

    //  // derivative of gradients at cell centers
    // dqf.malloc((mesh.Nelements+mesh.totalHaloPairs)*Nsfields*mesh.dim*mesh.Nfaces); 
    // o_dqf = platform.malloc<dfloat>(dqf); 

    // derivative of gradients at cell centers
    dq.malloc((mesh.Nelements+mesh.totalHaloPairs)*Nsfields*mesh.dim); 
    o_dq = platform.malloc<dfloat>(dq); 


  props["defines/" "s_Nvgeo"] = Nvgeo;
  props["defines/" "s_CXID"]  = CXID;
  props["defines/" "s_CYID"]  = CYID;
  props["defines/" "s_VID"]   = VID;
  
  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = STAB_DIR "/okl/";
  std::string oklFileSuffix = ".okl";
  std::string fileName, kernelName;

  fileName      = oklFilePrefix + "limiter" + oklFileSuffix;
  kernelName    = "projectCellAverage" + suffix;
  projectCellAverageKernel  = platform.buildKernel(fileName, kernelName, props);


  kernelName    = "getVertexValues";
  getVertexValuesKernel  = platform.buildKernel(fileName, kernelName, props);
}



void stab_t::stabApplyLimiter(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // detect elements
  Detect(o_Q, o_RHS, T); 

  // Project Solution to Cell Averages
  projectCellAverageKernel(mesh.Nelements, 
                           o_projectC0, 
                           o_Q,
                           o_qc,
                           o_qv); 

  ogs.GatherScatter(o_qv, Nsfields, ogs::Add, ogs::Sym);

  getVertexValuesKernel(mesh.Nelements, o_weight, o_qv); 

  // platform.linAlg().amx(mesh.Nelements*mesh.Nverts*Nsfields, 1.0, o_weight, o_qv); 

  // Exchange cell center data for reconstruction 
  mesh.halo.Exchange(o_qc, Nsfields); 
}




// Project to cell averages for limiting
void stab_t::limiterSetupOperators(){

  // Build gather scatter for vertex-based smmothing
  mesh_t meshC = mesh.SetupNewDegree(1);

  // Define wights for vertex smoothing
  bool verbose = false,  unique  = true; 
  ogs.Setup(mesh.Nelements*mesh.Nverts, meshC.globalIds, mesh.comm, 
                ogs::Signed, ogs::Auto, unique, verbose, platform); 
  
  //  
  weight.malloc(mesh.Nelements*mesh.Nverts, 1.0);
  ogs.GatherScatter(weight,1, ogs::Add, ogs::Sym); 

  // invert weights for gs for vertex nodes only
  for(int i=0; i < mesh.Nelements*mesh.Nverts; i++ ){ 
    // printf(" %d %d %.4f  %.4f \n", i/3, i%3, weight[i], 1.0/weight[i]);
    weight[i] = 1./weight[i];

  }
  

  // Move to device
  o_weight = platform.malloc<dfloat>(weight); 



  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      // Project to Cell Center Values
      projectC0.malloc(mesh.Np, 0.0);
      // // Compute Operator for Taylor Expansion x_mean - x
      // dAVE.malloc(mesh.Np*mesh.Np, 0.0); 
      for(int j=0; j<mesh.Np; j++){
        dfloat sum=0.0; 
        for(int i=0; i<mesh.Np; i++){
          sum += mesh.MM[i*mesh.Np + j]; 
      }
      projectC0[j] = 0.5*sum; 

      // printf(" %.4e \n ", projectC0[j]);

      // // dAVE = eye(Np) - projectC0
      // for(int i=0; i<mesh.Np; i++){
      //   dAVE[i*mesh.Np + j ] = (i==j) ? (1.0- projectC0[j]):-projectC0[j] ; 
      //   }
      }
      break; 
    case Mesh::TETRAHEDRA:
      projectC0.malloc(mesh.Np, 0.0);
      for(int j=0; j<mesh.Np; j++){
        dfloat sum=0.0; 
        for(int i=0; i<mesh.Np; i++){
          sum += mesh.MM[i*mesh.Np + j]; 
      }
      projectC0[j] = 0.5*sum; 
      }
      break; 
    case Mesh::QUADRILATERALS:{
      const int Np1 = mesh.N + 1; 
      memory<dfloat> r1, V1, M1; 
      mesh.Nodes1D(mesh.N, r1);
      mesh.Vandermonde1D(mesh.N, r1, V1); 
      mesh.MassMatrix1D(Np1, V1, M1); 

      projectC0.malloc(Np1, 0.0);
      for(int j=0; j<Np1; j++){
        dfloat sum=0.0; 
        for(int i=0; i<Np1; i++){
          sum += M1[i*Np1 + j]; 
      }
      projectC0[j] = 0.5*sum; 
      }
    }break; 
    case Mesh::HEXAHEDRA:
      break; 
  }
// Move to device
o_projectC0 = platform.malloc<dfloat>(projectC0); 
}



// Compute required geometric factors
void stab_t::limiterGeometricFactorsTri2D(){

Nvgeo = 3; 
CXID = 0; CYID = 1; VID = 2;  
// CXID = 0; CYID = 1; SAID = 2;  

vgeo.malloc((mesh.Nelements+mesh.totalHaloPairs)*Nvgeo); 
lvmapM.malloc(mesh.Nfaces*mesh.NfaceVertices*mesh.Nelements); 
lvmapP.malloc(mesh.Nfaces*mesh.NfaceVertices*mesh.Nelements); 


// First Compute cell Centers Coordinates
for(int e=0; e<mesh.Nelements; e++){
  dfloat xc0 = 0.0, yc0 = 0.0; 
  for(int n=0; n<mesh.Np; n++){
    xc0 += projectC0[n]*mesh.x[e*mesh.Np + n]; 
    yc0 += projectC0[n]*mesh.y[e*mesh.Np + n]; 
  }
  vgeo[e*Nvgeo + CXID] = xc0; 
  vgeo[e*Nvgeo + CYID] = yc0; 
  vgeo[e*Nvgeo + VID] = 2.0*mesh.vgeo[e*mesh.Nvgeo + mesh.JID]/3.0; 
  //
  // get vertex coordinates for each face pair
  // const dlong base = e*mesh.Nfaces*mesh.NfaceVertices; 
  // const dlong  idM = e*mesh.Nfaces*mesh.Nfp ;  
  // for(int face=0; face<mesh.Nfaces; face++){
  //   // fnode1 =  mesh.faceNodes[face*mesh.Nfp + 0];   
  //   // fnode2 =  mesh.faceNodes[(face+1)*mesh.Nfp - 1];
  //   lvmapM[base+ face*mesh.NfaceVertices + 0] =  mesh.vmapM[idM + face*mesh.Nfp + 0]; 
  //   lvmapM[base+ face*mesh.NfaceVertices + 1] =  mesh.vmapM[idM + (face+1)*mesh.Nfp - 1]; 
  //   lvmapP[base+ face*mesh.NfaceVertices + 0] =  mesh.vmapP[idM + face*mesh.Nfp + 0]; 
  //   lvmapP[base+ face*mesh.NfaceVertices + 1] =  mesh.vmapP[idM + (face+1)*mesh.Nfp - 1]; 
  // }
}

  
// Setup OGS for Vertex only Gather-Scattwer Operation
mesh.halo.Exchange(vgeo, Nvgeo);
o_vgeo = platform.malloc<dfloat>(vgeo); 
// o_lvmapM = platform.malloc<dlong>(lvmapM); 
// o_lvmapP = platform.malloc<dlong>(lvmapP); 



}


// Compute required geometric factors
void stab_t::limiterGeometricFactorsQuad2D(){




}


// Compute required geometric factors
void stab_t::limiterGeometricFactorsTet3D(){




}


// Compute required geometric factors
void stab_t::limiterGeometricFactorsHex3D(){




}









} //namespace libp
