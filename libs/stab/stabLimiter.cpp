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

namespace libp {

void stab_t::stabSetupLimiter(){
  /*

  // Setup Operators
  limSetupOperators(); 
  // Setup Geometry
   switch (mesh.elementType) {
      case Mesh::TRIANGLES:
        limGeometricFactorsTri2D();
        break;
      case Mesh::QUADRILATERALS:
        // limGeometricFactorsQuad2D(); 
        LIBP_FORCE_ABORT("Limiter-Quad is not implemented yet");
        break;
      case Mesh::TETRAHEDRA:
        // limGeometricFactorsTet3D(); 
        LIBP_FORCE_ABORT("Limiter-Quad is not implemented yet");
        break;
      case Mesh::HEXAHEDRA:
        // limGeometricFactorsHex3D();
        LIBP_FORCE_ABORT("Limiter-Hex is not implemented yet");
        break;
    }

    vertexNodes.malloc(mesh.Nverts); 
    vertexNodes.copyFrom(mesh.vertexNodes); 
    o_vertexNodes = platform.malloc<int>(vertexNodes);

    // Gather scatter vertex nodes
    qv.malloc(mesh.Nelements*mesh.Nfaces*mesh.NfaceVertices*sNfields); 
    o_qv = platform.malloc<dfloat>(qv); 

    // printf("%d %d\n", mesh.Nfaces, mesh.NfaceVertices); 


    qc.malloc((mesh.Nelements + mesh.totalHaloPairs)*sNfields); 
    o_qc = platform.malloc<dfloat>(qc); 


     // derivative of gradients at cell centers
    dqf.malloc((mesh.Nelements+mesh.totalHaloPairs)*sNfields*mesh.dim*mesh.Nfaces); 
    o_dqf = platform.malloc<dfloat>(dqf); 

    // derivative of gradients at cell centers
    dq.malloc((mesh.Nelements+mesh.totalHaloPairs)*sNfields*mesh.dim); 
    o_dq = platform.malloc<dfloat>(dq); 


    visc.malloc((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*dNfields, 0.0); 
     o_visc = platform.malloc<dfloat>(visc); 

     // Moved to solver AK
    // int blockMax = 256;
    // if (platform.device.mode() == "CUDA") blockMax = 512;

    // printf("right here\n");
    // props["defines/" "s_Nvgeo"]= Nvgeo;
    // props["defines/" "s_CXID"]= CXID;
    // props["defines/" "s_CYID"]= CYID;
    // props["defines/" "s_SAID"]= SAID;
    // props["defines/" "s_Nverts"]= int(mesh.Nverts);
    */

    
}


// Project to cell averages for limiting
void stab_t::limSetupOperators(){

  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      // Project to Cell Center Values
      projectC0.malloc(mesh.Np, 0.0);
      // Compute Operator for Taylor Expansion x_mean - x
      dAVE.malloc(mesh.Np*mesh.Np, 0.0); 
      for(int j=0; j<mesh.Np; j++){
        dfloat sum=0.0; 
        for(int i=0; i<mesh.Np; i++){
          sum += mesh.MM[i*mesh.Np + j]; 
      }
      projectC0[j] = 0.5*sum; 

      // dAVE = eye(Np) - projectC0
      for(int i=0; i<mesh.Np; i++){
        dAVE[i*mesh.Np + j ] = i==j ? (1.0- projectC0[j]):-projectC0[j] ; 
      }
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
void stab_t::limGeometricFactorsTri2D(){

Nvgeo = 3; 
CXID = 0; CYID = 1; SAID = 2;  

vgeo.malloc((mesh.Nelements+mesh.totalHaloPairs)*Nvgeo); 

DX.malloc(mesh.Nelements*mesh.Np*mesh.dim); 

// First Compute cell Centers Coordinates
for(int e=0; e<mesh.Nelements; e++){
  dfloat xc0 = 0.0, yc0 = 0.0; 
  for(int n=0; n<mesh.Np; n++){
    xc0 += projectC0[n]*mesh.x[e*mesh.Np + n]; 
    yc0 += projectC0[n]*mesh.y[e*mesh.Np + n]; 
  }
  vgeo[e*Nvgeo + CXID] = xc0; 
  vgeo[e*Nvgeo + CYID] = yc0; 
}

// Compute constructed areas
for(int e=0; e<mesh.Nelements; e++){
    vgeo[e*Nvgeo + SAID] = 2.0*mesh.vgeo[e*mesh.Nvgeo + mesh.JID]/3.0; 
}

for(int e=0; e<mesh.Nelements; e++){
 for(int n=0; n<mesh.Np; n++){
  dfloat sumx = 0. , sumy=0.; 
  for(int i=0; i<mesh.Np; i++){
    sumx += dAVE[n*mesh.Np+ i]*mesh.x[e*mesh.Np + i]; 
    sumy += dAVE[n*mesh.Np+ i]*mesh.y[e*mesh.Np + i]; 
    }
    DX[e*mesh.Np*mesh.dim + 0*mesh.Np + n] = sumx; 
    DX[e*mesh.Np*mesh.dim + 1*mesh.Np + n] = sumy; 
  }
}
o_DX = platform.malloc<dfloat>(DX); 

// Setup OGS for Vertex only Gather-Scattwer Operation
mesh.halo.Exchange(vgeo, Nvgeo);
o_vgeo = platform.malloc<dfloat>(vgeo); 
}


// Compute required geometric factors
void stab_t::limGeometricFactorsQuad2D(){




}


// Compute required geometric factors
void stab_t::limGeometricFactorsTet3D(){




}


// Compute required geometric factors
void stab_t::limGeometricFactorsHex3D(){




}









} //namespace libp
