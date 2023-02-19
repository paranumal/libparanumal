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
void SEMFEMPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  pfloat zero = 0, one = 1;
  
  linAlg_t& linAlg = elliptic.platform.linAlg();

  if (mesh.elementType==Mesh::TRIANGLES) {
    dlong Ncols = parAlmond.getNumCols(0);
    deviceMemory<pfloat> o_GrFEM = elliptic.platform.reserve<pfloat>(Ncols);
    deviceMemory<pfloat> o_GzFEM = elliptic.platform.reserve<pfloat>(Ncols);

    // Mr = invDegree.*r
    linAlg.amxpy(elliptic.Ndofs, one, elliptic.o_weightG, o_r, zero, o_Mr);

    deviceMemory<pfloat> o_rFEM = elliptic.platform.reserve<pfloat>(mesh.Nelements*mesh.NpFEM);
    elliptic.gHalo.Exchange(o_Mr, 1);
    SEMFEMInterpKernel(mesh.Nelements,
                       elliptic.o_GlobalToLocal,
                       mesh.o_pfloat_SEMFEMAnterp,
                       o_Mr, o_rFEM);
    FEMogs.Gather(o_GrFEM, o_rFEM, 1, ogs::Add, ogs::Trans);
    o_rFEM.free();

    parAlmond.Operator(o_GrFEM, o_GzFEM);

    deviceMemory<pfloat> o_MrL = elliptic.platform.reserve<pfloat>(mesh.Np*mesh.Nelements);
    FEMgHalo.Exchange(o_GzFEM, 1);
    SEMFEMAnterpKernel(mesh.Nelements,
                       o_FEMGlobalToLocal,
                       mesh.o_pfloat_SEMFEMAnterp,
                       o_GzFEM, o_MrL);
    elliptic.ogsMasked.Gather(o_Mr, o_MrL, 1, ogs::Add, ogs::Trans);
    o_MrL.free();

    // Mr = invDegree.*Mr
    linAlg.amx(elliptic.Ndofs, one, elliptic.o_weightG, o_Mr);

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

  /* Build halo exchange for gathered ordering */
  FEMgHalo.SetupFromGather(FEMogs);

  FEMGlobalToLocal.malloc(mesh.Nelements*mesh.NpFEM);
  FEMogs.SetupGlobalToLocalMapping(FEMGlobalToLocal);

  o_FEMGlobalToLocal = elliptic.platform.malloc<dlong>(FEMGlobalToLocal);

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

  memory<pfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);

  parAlmond.Report();

  if (mesh.elementType==Mesh::TRIANGLES) {
    // build interp and anterp
    memory<pfloat> pfloat_SEMFEMAnterp(mesh.NpFEM*mesh.Np);
    memory<pfloat> pfloat_SEMFEMInterp(mesh.NpFEM*mesh.Np);
    for(int n=0;n<mesh.NpFEM;++n){
      for(int m=0;m<mesh.Np;++m){
        pfloat_SEMFEMAnterp[n+m*mesh.NpFEM] = mesh.SEMFEMInterp[n*mesh.Np+m];
	pfloat_SEMFEMInterp[n+m*mesh.NpFEM] = mesh.SEMFEMInterp[n+m*mesh.NpFEM];
      }
    }

    mesh.o_pfloat_SEMFEMInterp = elliptic.platform.malloc<pfloat>(pfloat_SEMFEMInterp);
    mesh.o_pfloat_SEMFEMAnterp = elliptic.platform.malloc<pfloat>(pfloat_SEMFEMAnterp);

    //build kernels
    properties_t kernelInfo = mesh.props;

    kernelInfo["defines/" "p_Np"]= mesh.Np;
    kernelInfo["defines/" "p_NpFEM"]= mesh.NpFEM;

    int NblockV = std::max(256/mesh.NpFEM, 1);
    kernelInfo["defines/" "p_NblockV"]= NblockV;

    kernelInfo["defines/" "dfloat"]= pfloatString;

    SEMFEMInterpKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
                                     "ellipticSEMFEMInterp", kernelInfo);

    SEMFEMAnterpKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
                                     "ellipticSEMFEMAnterp", kernelInfo);
  }
}
