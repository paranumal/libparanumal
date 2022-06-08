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

#include "ellipticPrecon.hpp"

// Cast problem into spectrally-equivalent N=1 FEM space and precondition with AMG
void SEMFEMPrecon::Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr) {

  linAlg_t& linAlg = elliptic.platform.linAlg();

  if (mesh.elementType==Mesh::TRIANGLES) {

    // Mr = invDegree.*r
    linAlg.amxpy(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_r, 0.0, o_Mr);

    elliptic.ogsMasked.Scatter(o_MrL, o_Mr, 1, ogs::NoTrans);
    SEMFEMInterpKernel(mesh.Nelements, mesh.o_SEMFEMAnterp, o_MrL, o_rFEM);
    FEMogs.Gather(o_GrFEM, o_rFEM, 1, ogs::Add, ogs::Trans);

    parAlmond.Operator(o_GrFEM, o_GzFEM);

    FEMogs.Scatter(o_zFEM, o_GzFEM, 1, ogs::NoTrans);
    SEMFEMAnterpKernel(mesh.Nelements, mesh.o_SEMFEMAnterp, o_zFEM, o_MrL);
    elliptic.ogsMasked.Gather(o_Mr, o_MrL, 1, ogs::Add, ogs::Trans);

    // Mr = invDegree.*Mr
    linAlg.amx(elliptic.Ndofs, 1.0, elliptic.o_weightG, o_Mr);

  } else {
    //pass to parAlmond
    parAlmond.Operator(o_r, o_Mr);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

SEMFEMPrecon::SEMFEMPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings, mesh.comm) {

  //sanity checking
  LIBP_ABORT("SEMFEM is supported for CONTINUOUS only",
             !settings.compareSetting("DISCRETIZATION", "CONTINUOUS"));

  //make a low-order fem mesh from the sem mesh (also return globalIds of the enriched sem nodes, and faceNode mapping)
  memory<hlong> globalIds;
  memory<int> mapB;
  femMesh = mesh.SetupSEMFEM(globalIds, mapB);

  //use the BCs to make a maskedGlobalIds array
  dlong Ntotal = mesh.NpFEM*mesh.Nelements;
  memory<hlong> maskedGlobalIds(Ntotal);
  maskedGlobalIds.copyFrom(globalIds, Ntotal);

  if (mesh.elementType==Mesh::TRIANGLES) { //build a new mask for NpFEM>Np node sets
    //translate the node-wise bc flag
    for (int n=0;n<Ntotal;n++) {
      int bc = mapB[n];
      if (bc>0) {
        int BC = elliptic.BCType[bc];     //translate mesh's boundary flag
        mapB[n] = BC;  //record it

        if (mapB[n] == 1) maskedGlobalIds[n] = 0;   //Dirichlet boundary
      }
    }
  } else {
    //mask using the original mask
    for (dlong n=0;n<elliptic.Nmasked;n++)
      maskedGlobalIds[elliptic.maskIds[n]] = 0;
  }

  //build masked gs handle to gather from enriched sem nodes to assembled fem problem
  int verbose = 0;
  bool unique = true;
  FEMogs.Setup(Ntotal, maskedGlobalIds,
               mesh.comm, ogs::Signed, ogs::Auto,
               unique, verbose, elliptic.platform);

  //make a map from the fem mesh's nodes to the (enriched) sem nodes
  memory<dlong> localIds(femMesh.Nelements*femMesh.Nverts);
  for(dlong e=0;e<mesh.Nelements;++e){
    for (int n=0;n<mesh.NelFEM;n++) {
      dlong id[femMesh.Nverts];

      //local ids in the subelement fem grid
      for (int i=0;i<mesh.Nverts;i++)
        id[i] = e*mesh.NpFEM + mesh.FEMEToV[n*femMesh.Nverts+i];

      dlong femId = e*mesh.NelFEM*femMesh.Nverts+n*mesh.Nverts;
      switch(mesh.elementType){
      case Mesh::TRIANGLES:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[2];
        break;
      case Mesh::QUADRILATERALS:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[3];  //need to swap this as the Np nodes are ordered [0,1,3,2] in a degree 1 element
        localIds[femId+3] = id[2];
        break;
      case Mesh::TETRAHEDRA:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[2];
        localIds[femId+3] = id[3];
        break;
      case Mesh::HEXAHEDRA:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[3];  //need to swap this as the Np nodes are ordered [0,1,3,2,4,5,7,6] in a degree 1 element
        localIds[femId+3] = id[2];
        localIds[femId+4] = id[4];
        localIds[femId+5] = id[5];
        localIds[femId+6] = id[7];
        localIds[femId+7] = id[6];
        break;
      }
    }
  }

  //make a fem elliptic solver
  femElliptic.platform = elliptic.platform;
  femElliptic.mesh = femMesh;
  femElliptic.settings = elliptic.settings;
  femElliptic.lambda = elliptic.lambda;

  femElliptic.ogsMasked = FEMogs; //only for getting Ngather when building matrix

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = FEMogs.Ngather;

  // create a global numbering system
  memory<hlong> globalIds2(Ngather);

  // every gathered degree of freedom has its own global id
  hlong globalGatherOffset = Ngather;
  mesh.comm.Scan(Ngather, globalGatherOffset);
  globalGatherOffset = globalGatherOffset - Ngather;

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<FEMogs.Ngather;n++) {
    globalIds2[n] = n + globalGatherOffset;
  }

  //scatter this numbering to the original nodes
  Ntotal = mesh.NpFEM*mesh.Nelements;
  memory<hlong> maskedGlobalNumbering(Ntotal, -1);
  FEMogs.Scatter(maskedGlobalNumbering, globalIds2, 1, ogs::NoTrans);

  //transfer the consecutive global numbering to the fem mesh
  Ntotal = femMesh.Np*femMesh.Nelements;
  femElliptic.maskedGlobalNumbering.malloc(Ntotal);

  for (dlong e=0;e<femMesh.Nelements;e++) {
    for (int n=0;n<femMesh.Np;n++) {
      dlong id = e*femMesh.Np + n;
      dlong localId = localIds[id];
      femElliptic.maskedGlobalNumbering[id] = maskedGlobalNumbering[localId];
    }
  }

  //finally, build the fem matrix and pass to parAlmond
  if (mesh.rank==0){
    printf("-----------------------------Multigrid AMG Setup--------------------------------------------\n");
  }
  parAlmond::parCOO A(elliptic.platform, femMesh.comm);
  femElliptic.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = femMesh.rank;
  int size = femMesh.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = static_cast<dlong>(A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);

  memory<dfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);

  parAlmond.Report();

  if (mesh.elementType==Mesh::TRIANGLES) {
    // build interp and anterp
    memory<dfloat> SEMFEMAnterp(mesh.NpFEM*mesh.Np);
    for(int n=0;n<mesh.NpFEM;++n){
      for(int m=0;m<mesh.Np;++m){
        SEMFEMAnterp[n+m*mesh.NpFEM] = mesh.SEMFEMInterp[n*mesh.Np+m];
      }
    }

    mesh.o_SEMFEMInterp = elliptic.platform.malloc<dfloat>(mesh.SEMFEMInterp);
    mesh.o_SEMFEMAnterp = elliptic.platform.malloc<dfloat>(SEMFEMAnterp);

    memory<dfloat> dummy(mesh.Nelements*mesh.NpFEM,0.0); //need this to avoid uninitialized memory warnings
    o_rFEM = elliptic.platform.malloc<dfloat>(dummy);
    o_zFEM = elliptic.platform.malloc<dfloat>(dummy);

    dlong Ncols = parAlmond.getNumCols(0);
    dummy.malloc(Ncols,0.0);
    o_GrFEM = elliptic.platform.malloc<dfloat>(dummy);
    o_GzFEM = elliptic.platform.malloc<dfloat>(dummy);

    o_MrL = elliptic.platform.malloc<dfloat>(mesh.Np*mesh.Nelements);

    //build kernels
    properties_t kernelInfo = mesh.props;

    kernelInfo["defines/" "p_Np"]= mesh.Np;
    kernelInfo["defines/" "p_NpFEM"]= mesh.NpFEM;

    int NblockV = std::max(256/mesh.NpFEM, 1);
    kernelInfo["defines/" "p_NblockV"]= NblockV;

    SEMFEMInterpKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
                                     "ellipticSEMFEMInterp", kernelInfo);

    SEMFEMAnterpKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
                                     "ellipticSEMFEMAnterp", kernelInfo);
  }
}
