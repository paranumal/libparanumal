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

#include "mesh.hpp"

namespace libp {

mesh_t mesh_t::SetupSEMFEM(memory<hlong>& globalIds_,
                           memory<int>& mapB_){

  //partially assembled fem mesh (result of projecting sem element to larger space)
  mesh_t pmesh=*this;

  //setup the intermediate mesh for tris and tets
  if (elementType==Mesh::TRIANGLES) {
    /* SEMFEM data */
    SEMFEMNodesTri2D(N, NpFEM, rFEM, sFEM);
    SEMFEMEToVTri2D(N, NelFEM, FEMEToV);

    SEMFEMInterpMatrixTri2D(N, r, s, rFEM, sFEM, SEMFEMInterp);

    //set semfem nodes as the grid points
    pmesh.Np = NpFEM;
    pmesh.r  = rFEM;
    pmesh.s  = sFEM;

    //count number of face nodes in the semfem element
    dfloat NODETOL = 1e-6;
    pmesh.Nfp=0;
    for (int n=0;n<pmesh.Np;n++)
      if (std::abs(pmesh.s[n]+1)<NODETOL) pmesh.Nfp++;

    //remake the faceNodes array
    pmesh.faceNodes.malloc(Nfaces*pmesh.Nfp);
    int f0=0, f1=0, f2=0;
    for (int n=0;n<pmesh.Np;n++) {
      if (std::abs(pmesh.s[n]+1)<NODETOL)          pmesh.faceNodes[0*pmesh.Nfp+f0++] = n;
      if (std::abs(pmesh.r[n]+pmesh.s[n])<NODETOL) pmesh.faceNodes[1*pmesh.Nfp+f1++] = n;
      if (std::abs(pmesh.r[n]+1)<NODETOL)          pmesh.faceNodes[2*pmesh.Nfp+f2++] = n;
    }

    //remake vertexNodes array
    pmesh.vertexNodes.malloc(Nverts);
    for(int n=0;n<pmesh.Np;++n){
      if( (pmesh.r[n]+1)*(pmesh.r[n]+1)+(pmesh.s[n]+1)*(pmesh.s[n]+1)<NODETOL)
        pmesh.vertexNodes[0] = n;
      if( (pmesh.r[n]-1)*(pmesh.r[n]-1)+(pmesh.s[n]+1)*(pmesh.s[n]+1)<NODETOL)
        pmesh.vertexNodes[1] = n;
      if( (pmesh.r[n]+1)*(pmesh.r[n]+1)+(pmesh.s[n]-1)*(pmesh.s[n]-1)<NODETOL)
        pmesh.vertexNodes[2] = n;
    }

    // compute physical (x,y) locations FEM vertices
    pmesh.PhysicalNodes();

    // connect face nodes (find trace indices)
    pmesh.ConnectFaceNodes();

    // make a global indexing
    pmesh.ConnectNodes();
    //pmesh.globalIds is now populated
    //pmesh.mapB is now populated
  } else if (elementType==Mesh::QUADRILATERALS) {
    NpFEM = Np;
    NelFEM = N*N;
    SEMFEMEToVQuad2D(N, FEMEToV);
  } else if (elementType==Mesh::TETRAHEDRA){
    NpFEM = Np;
    NelFEM = N*N*N;
    SEMFEMEToVTet3D(N, FEMEToV);
  } else { //Mesh::HEXAHEDRA
    NpFEM = Np;
    NelFEM = N*N*N;
    SEMFEMEToVHex3D(N, FEMEToV);
  }

  //need to return this data
  globalIds_ = pmesh.globalIds;
  mapB_ = pmesh.mapB;

  //now build the full degree 1 fem mesh
  mesh_t femMesh=*this;

  femMesh.N = 1; //degree of fem approximation

  /* allocate space for node coordinates */
  femMesh.Nelements = NelFEM*Nelements;
  dlong NFEMverts = femMesh.Nelements*Nverts;
  femMesh.EToV.malloc(NFEMverts);
  femMesh.EX.malloc(NFEMverts);
  femMesh.EY.malloc(NFEMverts);
  if (dim==3)
    femMesh.EZ.malloc(NFEMverts);

  for(dlong e=0;e<Nelements;++e){
    for (int n=0;n<NelFEM;n++) {
      dlong femId = e*NelFEM*Nverts+n*Nverts;

      for (int i=0;i<Nverts;i++) {
        //local ids in the subelement fem grid
        dlong id = e*NpFEM + FEMEToV[n*Nverts+i];

        /* read vertex triplet for triangle */
        femMesh.EToV[femId+i] = pmesh.globalIds[id];

        femMesh.EX[femId+i] = pmesh.x[id];
        femMesh.EY[femId+i] = pmesh.y[id];
        if (dim==3)
          femMesh.EZ[femId+i] = pmesh.z[id];
      }
    }
  }

  // load reference (r,s) element nodes
  femMesh.ReferenceNodes();

  // connect elements using parallel sort
  femMesh.Connect();

  //identify the nodes on the SEMFEM element faces
  memory<int> faceFlag(pmesh.Np*Nfaces, 0);
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<pmesh.Nfp;n++) {
      int id = pmesh.faceNodes[f*pmesh.Nfp+n];
      faceFlag[f*pmesh.Np + id] = 1; //flag the nodes on this face
    }
  }

  //map from faces of fem sub-elements to the macro element face number
  memory<int> femFaceMap(NelFEM*femMesh.Nfaces, 0);
  for (int n=0;n<NelFEM*femMesh.Nfaces;n++) femFaceMap[n] = -1;

  for (int n=0;n<NelFEM;n++) {
    for (int f=0;f<femMesh.Nfaces;f++) {

      for (int face=0; face<Nfaces;face++) {

        //count the nodes on this face which are on a macro face
        int NvertsOnFace = 0;
        for (int i=0;i<femMesh.Nfp;i++){
          int id = femMesh.faceNodes[f*femMesh.Nfp+i];
          int v  = FEMEToV[n*Nverts+id];
          NvertsOnFace += faceFlag[face*pmesh.Np + v];
        }
        if (NvertsOnFace == femMesh.Nfp)
          femFaceMap[n*femMesh.Nfaces+f] = face; //on macro face
      }
    }
  }

  //fill the boundary flag array from the original EToB
  femMesh.EToB.malloc(femMesh.Nelements*femMesh.Nfaces, 0);
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<NelFEM;n++) {
      for (int f=0;f<femMesh.Nfaces;f++) {
        int face = femFaceMap[n*femMesh.Nfaces+f];
        if (face>-1) {
          femMesh.EToB[(e*NelFEM +n)*femMesh.Nfaces +f] = EToB[e*Nfaces + face];
        }
      }
    }
  }

  // set up halo exchange info for MPI (do before connect face nodes)
  femMesh.HaloSetup();

  // connect face vertices
  femMesh.ConnectFaceVertices();

  // connect face nodes (find trace indices)
  femMesh.ConnectFaceNodes();

  // make global indexing
  femMesh.ConnectNodes();

  // compute physical (x,y) locations of the element nodes
  femMesh.PhysicalNodes();

  // compute geometric factors
  femMesh.GeometricFactors();

  // compute surface geofacs
  // femMesh.SurfaceGeometricFactors();

  // label local/global gather elements
  femMesh.GatherScatterSetup();

  return femMesh;
}

} //namespace libp
