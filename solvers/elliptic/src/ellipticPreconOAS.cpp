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

#include "ellipticPrecon.hpp"

// Overlapping additive Schwarz with patch problems consisting of the
//  entire local mesh + 1 ring overlap, solved with a local multigrid
//  precon and coarse problem consisting of the global degree 1
//  problem, solved with parAlmond
void OASPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  if (mesh.N>1) {
    if (elliptic.disc_c0) {
      //Scatter to localDof ordering, exchange ring halo,
      // then compress the ring mesh to globalDofs order.
      // TODO: Technically, these steps couple be fused
      // to a single operation, but currently theres no easy way
      // as the ordering of globalDofs between the original mesh
      // partition and the ring mesh could be different
      elliptic.ogsMasked->Scatter(o_rPatchL, o_r, ogs_dfloat, ogs_add, ogs_notrans);
      mesh.ringHalo->Exchange(o_rPatchL, mesh.Np, ogs_dfloat);
      ellipticPatch->ogsMasked->Gather(o_rPatch, o_rPatchL, ogs_dfloat, ogs_add, ogs_notrans);
    } else {
      o_rPatch.copyFrom(o_r, elliptic.Ndofs*sizeof(dfloat));
      mesh.ringHalo->Exchange(o_rPatch, mesh.Np, ogs_dfloat);
    }

    //Apply local patch precon
    preconPatch->Operator(o_rPatch, o_zPatch);

    //Coarsen problem to N=1 and pass to parAlmond
    // TODO: This is blocking due to H<->D transfers.
    //       Should modify precons so size=1 is non-blocking
    level->coarsen(o_r, o_rC);

    parAlmond.Operator(o_rC, o_zC);

    //Add contributions from all patches together
    if (elliptic.disc_c0) {
      dlong Ntotal=mesh.Nelements*mesh.Np;

      ellipticPatch->ogsMasked->Scatter(o_zPatchL, o_zPatch, ogs_dfloat, ogs_add, ogs_notrans);
      ogsMaskedRing->GatherScatter(o_zPatchL, ogs_dfloat, ogs_add, ogs_sym);

      // Weight by overlap degree, zPatch = patchWeight*zPatch
      elliptic.linAlg.amx(Ntotal, 1.0, o_patchWeight, o_zPatchL);

      elliptic.ogsMasked->Gather(o_Mr, o_zPatchL, ogs_dfloat, ogs_add, ogs_notrans);

    } else {
      mesh.ringHalo->Combine(o_zPatch, mesh.Np, ogs_dfloat);

      // Weight by overlap degree, Mr = patchWeight*zPatch
      elliptic.linAlg.amxpy(elliptic.Ndofs, 1.0, o_patchWeight, o_zPatch, 0.0, o_Mr);
    }

    // Add prologatated coarse solution
    level->prolongate(o_zC, o_Mr);
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
    meshPatch = mesh.SetupRingPatch();
    ellipticPatch = elliptic.SetupRingPatch(*meshPatch);
    preconPatch = new MultiGridPrecon(*ellipticPatch);

    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      rPatchL = (dfloat*) calloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),sizeof(dfloat));
      zPatchL = (dfloat*) calloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),sizeof(dfloat));

      o_rPatchL = elliptic.platform.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), rPatchL);
      o_zPatchL = elliptic.platform.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), zPatchL);
    }

    rPatch  = (dfloat*) calloc(ellipticPatch->Ndofs,sizeof(dfloat));
    zPatch  = (dfloat*) calloc(ellipticPatch->Ndofs,sizeof(dfloat));
    o_rPatch = elliptic.platform.malloc(ellipticPatch->Ndofs*sizeof(dfloat), rPatch);
    o_zPatch = elliptic.platform.malloc(ellipticPatch->Ndofs*sizeof(dfloat), zPatch);

    //compute patch overlap weighting
    patchWeight = (dfloat*) malloc(meshPatch->Nelements*meshPatch->Np*sizeof(dfloat));
    for (int i=0;i<meshPatch->Nelements*meshPatch->Np;i++)
      patchWeight[i] = 1.0;

    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      //share the masked version of the global id numbering
      hlong *maskedRingGlobalIds = (hlong *) calloc(meshPatch->Nelements*meshPatch->Np,sizeof(hlong));
      memcpy(maskedRingGlobalIds, elliptic.maskedGlobalIds, mesh.Nelements*mesh.Np*sizeof(hlong));
      mesh.ringHalo->Exchange(maskedRingGlobalIds, mesh.Np, ogs_hlong);

      //mask ring
      for (dlong n=0;n<ellipticPatch->Nmasked;n++)
        maskedRingGlobalIds[ellipticPatch->maskIds[n]] = 0;

      //use the masked ids to make another gs handle
      int verbose = 0;
      ogsMaskedRing = ogs_t::Setup(meshPatch->Nelements*meshPatch->Np, maskedRingGlobalIds,
                                   mesh.comm, verbose, elliptic.platform);
      free(maskedRingGlobalIds);

      //determine overlap of each node with masked ogs
      ogsMaskedRing->GatherScatter(patchWeight, ogs_dfloat, ogs_add, ogs_sym);

    } else {
      //determine overlap by combining halos
      mesh.ringHalo->Combine(patchWeight, mesh.Np, ogs_dfloat);
    }

    //invert
    for (int i=0;i<meshPatch->Nelements*meshPatch->Np;i++)
      patchWeight[i] = (patchWeight[i] > 0.0) ? 1.0/patchWeight[i] : 0.0;

    o_patchWeight = elliptic.platform.malloc(meshPatch->Nelements*meshPatch->Np*sizeof(dfloat), patchWeight);
  }

  //build the coarse precon
  int Nc = 1;  //hard code
  int NpCoarse = mesh.Np;
  switch(mesh.elementType){
    case TRIANGLES:
      NpCoarse = ((Nc+1)*(Nc+2))/2; break;
    case QUADRILATERALS:
      NpCoarse = (Nc+1)*(Nc+1); break;
    case TETRAHEDRA:
      NpCoarse = ((Nc+1)*(Nc+2)*(Nc+3))/6; break;
    case HEXAHEDRA:
      NpCoarse = (Nc+1)*(Nc+1)*(Nc+1); break;
  }

  //build mesh and elliptic objects for this degree
  mesh_t &meshC = mesh.SetupNewDegree(Nc);
  elliptic_t &ellipticC = elliptic.SetupNewDegree(meshC);

  //build full A matrix and pass to parAlmond
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
  dfloat *null = (dfloat *) malloc(numLocalRows*sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1.0/sqrt(TotalRows);

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, ellipticC.allNeumann, null,ellipticC.allNeumannPenalty);

  if (mesh.N>1) {
    //make an MG level to get prologation and coarsener
    dlong Nrows, Ncols;
    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      Nrows = elliptic.ogsMasked->Ngather;
      Ncols = Nrows + elliptic.ogsMasked->NgatherHalo;
    } else {
      Nrows = mesh.Nelements*mesh.Np;
      Ncols = Nrows + mesh.totalHaloPairs*mesh.Np;
    }

    level = new MGLevel(elliptic, Nrows, Ncols, Nc, NpCoarse);
    level->meshC = &meshC;
    level->ogsMaskedC = ellipticC.ogsMasked;

    //coarse buffers
    Ncols = parAlmond.getNumCols(0);
    rC = (dfloat*) calloc(Ncols,sizeof(dfloat));
    zC = (dfloat*) calloc(Ncols,sizeof(dfloat));
    o_rC = elliptic.platform.malloc(Ncols*sizeof(dfloat), rC);
    o_zC = elliptic.platform.malloc(Ncols*sizeof(dfloat), zC);
  }

  //report
  parAlmond.Report();
}

OASPrecon::~OASPrecon() {
  if (mesh.N>1) {
    delete preconPatch;
    if (mesh.size>1) delete ellipticPatch;
    if (mesh.size>1) delete meshPatch;

    delete &(level->elliptic);
    if (level->mesh.ogs) level->mesh.ogs->Free();
    delete level;
  }
}