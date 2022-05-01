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

// Overlapping additive Schwarz with patch problems consisting of the
//  entire local mesh + 1 ring overlap, solved with a local multigrid
//  precon and coarse problem consisting of the global degree 1
//  problem, solved with parAlmond
void OASPrecon::Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr) {

  if (mesh.N>1) {
    if (elliptic.disc_c0) {
      //Scatter to localDof ordering, exchange ring halo,
      // then compress the ring mesh to globalDofs order.
      // TODO: Technically, these steps couple be fused
      // to a single operation, but currently theres no easy way
      // as the ordering of globalDofs between the original mesh
      // partition and the ring mesh could be different
      elliptic.ogsMasked.Scatter(o_rPatchL, o_r, 1, ogs::NoTrans);
      mesh.ringHalo.Exchange(o_rPatchL, mesh.Np);
      ellipticPatch.ogsMasked.Gather(o_rPatch, o_rPatchL, 1, ogs::Add, ogs::NoTrans);
    } else {
      o_rPatch.copyFrom(o_r, elliptic.Ndofs);
      mesh.ringHalo.Exchange(o_rPatch, mesh.Np);
    }

    //Apply local patch precon
    preconPatch.Operator(o_rPatch, o_zPatch);

    //Coarsen problem to N=1 and pass to parAlmond
    // TODO: This is blocking due to H<->D transfers.
    //       Should modify precons so size=1 is non-blocking
    level.coarsen(o_r, o_rC);

    parAlmond.Operator(o_rC, o_zC);

    linAlg_t& linAlg = elliptic.platform.linAlg();

    //Add contributions from all patches together
    if (elliptic.disc_c0) {
      dlong Ntotal=mesh.Nelements*mesh.Np;

      ellipticPatch.ogsMasked.Scatter(o_zPatchL, o_zPatch, 1, ogs::NoTrans);
      ogsMaskedRing.GatherScatter(o_zPatchL, 1, ogs::Add, ogs::Sym);

      // Weight by overlap degree, zPatch = patchWeight*zPatch
      linAlg.amx(Ntotal, 1.0, o_patchWeight, o_zPatchL);

      elliptic.ogsMasked.Gather(o_Mr, o_zPatchL, 1, ogs::Add, ogs::NoTrans);

    } else {
      mesh.ringHalo.Combine(o_zPatch, mesh.Np);

      // Weight by overlap degree, Mr = patchWeight*zPatch
      linAlg.amxpy(elliptic.Ndofs, 1.0, o_patchWeight, o_zPatch, 0.0, o_Mr);
    }

    // Add prologatated coarse solution
    level.prolongate(o_zC, o_Mr);
  } else {
    //if N=1 just call the coarse solver
    parAlmond.Operator(o_r, o_Mr);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

OASPrecon::OASPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings, mesh.comm) {

  //build the one ring mesh
  if (mesh.N>1) {
    if (Comm::World().rank()==0){
      printf("-----------------------------Multigrid Degree %2d Patch--------------------------------------\n", mesh.N);
    }
    meshPatch = mesh.SetupRingPatch();
    ellipticPatch = elliptic.SetupRingPatch(meshPatch);
    preconPatch.Setup<MultiGridPrecon>(ellipticPatch);

    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      rPatchL.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),0.0);
      zPatchL.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),0.0);

      o_rPatchL = elliptic.platform.malloc<dfloat>(rPatchL);
      o_zPatchL = elliptic.platform.malloc<dfloat>(zPatchL);
    }

    rPatch.malloc(ellipticPatch.Ndofs,0.0);
    zPatch.malloc(ellipticPatch.Ndofs,0.0);
    o_rPatch = elliptic.platform.malloc<dfloat>(rPatch);
    o_zPatch = elliptic.platform.malloc<dfloat>(zPatch);

    //compute patch overlap weighting
    patchWeight.malloc(meshPatch.Nelements*meshPatch.Np);
    for (int i=0;i<meshPatch.Nelements*meshPatch.Np;i++)
      patchWeight[i] = 1.0;

    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      //share the masked version of the global id numbering
      memory<hlong> maskedRingGlobalIds(meshPatch.Nelements*meshPatch.Np);
      maskedRingGlobalIds.copyFrom(elliptic.maskedGlobalIds, mesh.Nelements*mesh.Np);
      mesh.ringHalo.Exchange(maskedRingGlobalIds, mesh.Np);

      //mask ring
      for (dlong n=0;n<ellipticPatch.Nmasked;n++)
        maskedRingGlobalIds[ellipticPatch.maskIds[n]] = 0;

      //use the masked ids to make another gs handle
      int verbose = 0;
      bool unique = true; //flag a unique node in every gather node
      ogsMaskedRing.Setup(meshPatch.Nelements*meshPatch.Np,
                          maskedRingGlobalIds, mesh.comm,
                          ogs::Signed, ogs::Auto,
                          unique, verbose, elliptic.platform);

      //determine overlap of each node with masked ogs
      ogsMaskedRing.GatherScatter(patchWeight, 1, ogs::Add, ogs::Sym);

    } else {
      //determine overlap by combining halos
      mesh.ringHalo.Combine(patchWeight, mesh.Np);
    }

    //invert
    for (int i=0;i<meshPatch.Nelements*meshPatch.Np;i++)
      patchWeight[i] = (patchWeight[i] > 0.0) ? 1.0/patchWeight[i] : 0.0;

    o_patchWeight = elliptic.platform.malloc<dfloat>(patchWeight);
  }

  //build the coarse precon
  int Nc = 1;  //hard code
  int NpCoarse = mesh.Np;
  switch(mesh.elementType){
    case Mesh::TRIANGLES:
      NpCoarse = ((Nc+1)*(Nc+2))/2; break;
    case Mesh::QUADRILATERALS:
      NpCoarse = (Nc+1)*(Nc+1); break;
    case Mesh::TETRAHEDRA:
      NpCoarse = ((Nc+1)*(Nc+2)*(Nc+3))/6; break;
    case Mesh::HEXAHEDRA:
      NpCoarse = (Nc+1)*(Nc+1)*(Nc+1); break;
  }

  //build mesh and elliptic objects for this degree
  mesh_t meshC = mesh.SetupNewDegree(Nc);
  elliptic_t ellipticC = elliptic.SetupNewDegree(meshC);

  //build full A matrix and pass to parAlmond
  if (Comm::World().rank()==0){
    printf("-----------------------------Multigrid AMG Setup--------------------------------------------\n");
  }
  parAlmond::parCOO A(elliptic.platform, meshC.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    ellipticC.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    ellipticC.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = meshC.rank;
  int size = meshC.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = (dlong) (A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);
  memory<dfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, ellipticC.allNeumann, null,ellipticC.allNeumannPenalty);

  if (mesh.N>1) {
    //make an MG level to get prologation and coarsener
    dlong Nrows, Ncols;
    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      Nrows = elliptic.ogsMasked.Ngather;
      Ncols = Nrows + elliptic.gHalo.Nhalo;
    } else {
      Nrows = mesh.Nelements*mesh.Np;
      Ncols = Nrows + mesh.totalHaloPairs*mesh.Np;
    }

    level = MGLevel(elliptic, Nrows, Ncols, Nc, NpCoarse);
    level.meshC = meshC;
    level.ellipticC = ellipticC;

    //coarse buffers
    Ncols = parAlmond.getNumCols(0);
    rC.malloc(Ncols,0.0);
    zC.malloc(Ncols,0.0);
    o_rC = elliptic.platform.malloc<dfloat>(rC);
    o_zC = elliptic.platform.malloc<dfloat>(zC);
  }

  //report
  parAlmond.Report();
}
