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

#include "mesh.hpp"
#include "mesh2D.hpp"
#include "mesh3D.hpp"

mesh_t* mesh_t::SetupSEMFEM(hlong **globalIds_, int *Nfp_, int **faceNodes_){

  //partially assembled fem mesh (result of projecting sem element to larger space)
  mesh_t *pmesh=NULL;
  switch(elementType){
  //quads and hexes reuse the SEM ndoes for the FEM problem
  case QUADRILATERALS:
    pmesh = this; break;
  case TETRAHEDRA:
    pmesh = this; break;
  case HEXAHEDRA:
    pmesh = this; break;
  case TRIANGLES:
    if(dim==2)
      pmesh = new meshTri2D(device, comm, settings, props);
    else
      pmesh = new meshTri3D(device, comm, settings, props);
    break;
  }

  //setup the intermediate mesh for tris and tets
  if (elementType==TRIANGLES) {
    pmesh->dim           = dim;
    pmesh->elementType   = elementType;
    pmesh->Nverts        = Nverts;
    pmesh->Nfaces        = Nfaces;
    pmesh->NfaceVertices = NfaceVertices;
    pmesh->faceVertices  = faceVertices;

    /* SEMFEM data */
    SEMFEMNodesTri2D(N, &NpFEM, &rFEM, &sFEM);
    SEMFEMEToVTri2D(N, &NelFEM, &FEMEToV);

    SEMFEMInterp = (dfloat*) calloc(NpFEM*Np, sizeof(dfloat));
    SEMFEMInterpMatrixTri2D(N, Np, r, s, NpFEM, rFEM, sFEM, SEMFEMInterp);

    //set semfem nodes as the grid points
    pmesh->Np = NpFEM;
    pmesh->r  = rFEM;
    pmesh->s  = sFEM;

    //count number of face nodes in the semfem element
    dfloat NODETOL = 1e-6;
    pmesh->Nfp=0;
    for (int n=0;n<pmesh->Np;n++)
      if (fabs(pmesh->s[n]+1)<NODETOL) pmesh->Nfp++;

    //remake the faceNodes array
    pmesh->faceNodes = (int *) calloc(Nfaces*pmesh->Nfp,sizeof(int));
    int f0=0, f1=0, f2=0;
    for (int n=0;n<pmesh->Np;n++) {
      if (fabs(pmesh->s[n]+1)<NODETOL)           pmesh->faceNodes[0*pmesh->Nfp+f0++] = n;
      if (fabs(pmesh->r[n]+pmesh->s[n])<NODETOL) pmesh->faceNodes[1*pmesh->Nfp+f1++] = n;
      if (fabs(pmesh->r[n]+1)<NODETOL)           pmesh->faceNodes[2*pmesh->Nfp+f2++] = n;
    }

    //remake vertexNodes array
    pmesh->vertexNodes = (int*) calloc(Nverts, sizeof(int));
    for(int n=0;n<pmesh->Np;++n){
      if( (pmesh->r[n]+1)*(pmesh->r[n]+1)+(pmesh->s[n]+1)*(pmesh->s[n]+1)<NODETOL)
        pmesh->vertexNodes[0] = n;
      if( (pmesh->r[n]-1)*(pmesh->r[n]-1)+(pmesh->s[n]+1)*(pmesh->s[n]+1)<NODETOL)
        pmesh->vertexNodes[1] = n;
      if( (pmesh->r[n]+1)*(pmesh->r[n]+1)+(pmesh->s[n]-1)*(pmesh->s[n]-1)<NODETOL)
        pmesh->vertexNodes[2] = n;
    }

    // use existing mesh connectivity
    pmesh->Nnodes = Nnodes;
    pmesh->EX = EX; // coordinates of vertices for each element
    pmesh->EY = EY;
    pmesh->EZ = EZ;

    pmesh->Nelements = Nelements;
    pmesh->NelementsGlobal = NelementsGlobal;
    pmesh->EToV = EToV; // element-to-vertex connectivity
    pmesh->EToE = EToE; // element-to-element connectivity
    pmesh->EToF = EToF; // element-to-(local)face connectivity
    pmesh->EToP = EToP; // element-to-partition/process connectivity
    pmesh->EToB = EToB; // element-to-boundary condition type

    pmesh->elementInfo = elementInfo;

    pmesh->NboundaryFaces = NboundaryFaces;
    pmesh->boundaryInfo   = boundaryInfo;

    //use existing halo
    pmesh->halo = halo;
    pmesh->NinternalElements = NinternalElements;
    pmesh->NhaloElements = NhaloElements;
    pmesh->totalHaloPairs = totalHaloPairs;
    pmesh->internalElementIds = internalElementIds;
    pmesh->haloElementIds = haloElementIds;

    // compute physical (x,y) locations FEM vertices
    pmesh->PhysicalNodes();

    // connect face nodes (find trace indices)
    pmesh->ConnectFaceNodes();

    // make a global indexing
    pmesh->ParallelConnectNodes();
    //pmesh->globalIds is now populated
  }

  //need to return this data
  *globalIds_ = pmesh->globalIds;
  *Nfp_ = pmesh->Nfp;
  *faceNodes_ = pmesh->faceNodes;

  //now build the full degree 1 fem mesh
  mesh_t *femMesh=NULL;
  switch(elementType){
  case TRIANGLES:
    if(dim==2)
      femMesh = new meshTri2D(device, comm, settings, props);
    else
      femMesh = new meshTri3D(device, comm, settings, props);
    break;
  case QUADRILATERALS:
    if(dim==2)
      femMesh = new meshQuad2D(device, comm, settings, props);
    else
      femMesh = new meshQuad3D(device, comm, settings, props);
    NpFEM = Np;
    NelFEM = N*N;
    FEMEToV = (int*) malloc(NelFEM*Nverts*sizeof(int));
    SEMFEMEToVQuad2D(N, FEMEToV);
    break;
  case TETRAHEDRA:
    femMesh = new meshTet3D(device, comm, settings, props);
    NpFEM = Np;
    NelFEM = N*N*N;
    FEMEToV = (int*) malloc(NelFEM*Nverts*sizeof(int));
    SEMFEMEToVTet3D(N, FEMEToV);
    break;
  case HEXAHEDRA:
    femMesh = new meshHex3D(device, comm, settings, props);
    NpFEM = Np;
    NelFEM = N*N*N;
    FEMEToV = (int*) malloc(NelFEM*Nverts*sizeof(int));
    SEMFEMEToVHex3D(N, FEMEToV);
    break;
  }

  int femN = 1; //degree of fem approximation
  femMesh->dim           = dim;
  femMesh->elementType   = elementType;
  femMesh->Nverts        = Nverts;
  femMesh->Nfaces        = Nfaces;
  femMesh->NfaceVertices = NfaceVertices;
  femMesh->faceVertices  = faceVertices;

  /* allocate space for node coordinates */
  femMesh->Nelements = NelFEM*Nelements;
  dlong NFEMverts = femMesh->Nelements*Nverts;
  femMesh->EToV = (hlong*) calloc(NFEMverts, sizeof(hlong));
  femMesh->EX = (dfloat*) calloc(NFEMverts, sizeof(dfloat));
  femMesh->EY = (dfloat*) calloc(NFEMverts, sizeof(dfloat));
  if (dim==3)
    femMesh->EZ = (dfloat*) calloc(NFEMverts, sizeof(dfloat));

  for(dlong e=0;e<Nelements;++e){
    for (int n=0;n<NelFEM;n++) {
      dlong femId = e*NelFEM*Nverts+n*Nverts;

      for (int i=0;i<Nverts;i++) {
        //local ids in the subelement fem grid
        dlong id = e*NpFEM + FEMEToV[n*Nverts+i];

        /* read vertex triplet for triangle */
        femMesh->EToV[femId+i] = pmesh->globalIds[id];

        femMesh->EX[femId+i] = pmesh->x[id];
        femMesh->EY[femId+i] = pmesh->y[id];
        if (dim==3)
          femMesh->EZ[femId+i] = pmesh->z[id];
      }
    }
  }

  // connect elements using parallel sort
  femMesh->ParallelConnect();

  // load reference (r,s) element nodes
  femMesh->ReferenceNodes(femN);

  //identify the nodes on the SEMFEM element faces
  int *faceFlag = (int*) calloc(pmesh->Np*Nfaces,sizeof(int));
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<pmesh->Nfp;n++) {
      int id = pmesh->faceNodes[f*pmesh->Nfp+n];
      faceFlag[f*pmesh->Np + id] = 1; //flag the nodes on this face
    }
  }

  //map from faces of fem sub-elements to the macro element face number
  int *femFaceMap = (int*) calloc(NelFEM*femMesh->Nfaces,sizeof(int));
  for (int n=0;n<NelFEM*femMesh->Nfaces;n++) femFaceMap[n] = -1;

  for (int n=0;n<NelFEM;n++) {
    for (int f=0;f<femMesh->Nfaces;f++) {

      for (int face=0; face<Nfaces;face++) {

        //count the nodes on this face which are on a macro face
        int NvertsOnFace = 0;
        for (int i=0;i<femMesh->Nfp;i++){
          int id = femMesh->faceNodes[f*femMesh->Nfp+i];
          int v  = FEMEToV[n*Nverts+id];
          NvertsOnFace += faceFlag[face*pmesh->Np + v];
        }
        if (NvertsOnFace == femMesh->Nfp)
          femFaceMap[n*femMesh->Nfaces+f] = face; //on macro face
      }
    }
  }

  //fill the boundary flag array from the original EToB
  femMesh->EToB = (int*) calloc(femMesh->Nelements*femMesh->Nfaces, sizeof(int));
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<NelFEM;n++) {
      for (int f=0;f<femMesh->Nfaces;f++) {
        int face = femFaceMap[n*femMesh->Nfaces+f];
        if (face>-1) {
          femMesh->EToB[(e*NelFEM +n)*femMesh->Nfaces +f] = EToB[e*Nfaces + face];
        }
      }
    }
  }
  free(faceFlag);
  free(femFaceMap);

  // set up halo exchange info for MPI (do before connect face nodes)
  femMesh->HaloSetup();

  // compute physical (x,y) locations of the element nodes
  femMesh->PhysicalNodes();

  // compute geometric factors
  femMesh->GeometricFactors();

  // connect face nodes (find trace indices)
  // femMesh->ConnectFaceNodes();

  // compute surface geofacs
  // femMesh->SurfaceGeometricFactors();

  // make a global indexing
  //femMesh->ParallelConnectNodes();

  // make an ogs operator and label local/global gather elements
  //femMesh->ParallelGatherScatterSetup();

  //dont need to setup occa buffers for this mesh
  // femMesh->OccaSetup();

  //just build the stiffness ops
  if (elementType==TRIANGLES) {
    //build stiffness matrices
    femMesh->Srr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Srs = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Ssr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Sss = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    for (int n=0;n<femMesh->Np;n++) {
      for (int m=0;m<femMesh->Np;m++) {
        for (int k=0;k<femMesh->Np;k++) {
          for (int l=0;l<femMesh->Np;l++) {
            femMesh->Srr[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
            femMesh->Srs[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
            femMesh->Ssr[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
            femMesh->Sss[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
          }
        }
      }
    }
  } else if (elementType==TETRAHEDRA) {
    //build stiffness matrices
    femMesh->Srr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Srs = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Srt = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Ssr = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Sss = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Sst = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Str = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Sts = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    femMesh->Stt = (dfloat *) calloc(femMesh->Np*femMesh->Np,sizeof(dfloat));
    for (int n=0;n<femMesh->Np;n++) {
      for (int m=0;m<femMesh->Np;m++) {
        for (int k=0;k<femMesh->Np;k++) {
          for (int l=0;l<femMesh->Np;l++) {
            femMesh->Srr[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
            femMesh->Srs[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
            femMesh->Srt[m+n*femMesh->Np] += femMesh->Dr[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dt[m+k*femMesh->Np];
            femMesh->Ssr[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
            femMesh->Sss[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
            femMesh->Sst[m+n*femMesh->Np] += femMesh->Ds[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dt[m+k*femMesh->Np];
            femMesh->Str[m+n*femMesh->Np] += femMesh->Dt[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dr[m+k*femMesh->Np];
            femMesh->Sts[m+n*femMesh->Np] += femMesh->Dt[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Ds[m+k*femMesh->Np];
            femMesh->Stt[m+n*femMesh->Np] += femMesh->Dt[n+l*femMesh->Np]*femMesh->MM[k+l*femMesh->Np]*femMesh->Dt[m+k*femMesh->Np];
          }
        }
      }
    }
  }

  return femMesh;
}
