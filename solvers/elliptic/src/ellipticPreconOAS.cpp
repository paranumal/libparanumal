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

  dlong Ntotal = mesh.Np*mesh.Nelements;

  if (mesh.N>1) {
    //Copy to patch buffer and queue to the ring exchange
    o_rPatch.copyFrom(o_r, Ntotal*sizeof(dfloat));
    mesh.ringHalo->Exchange(o_rPatch, mesh.Np, ogs_dfloat);

    //Apply local patch precon
    preconPatch->Operator(o_rPatch, o_zPatch);

    //Coarsen problem to N=1 and pass to parAlmond
    // TODO: This is blocking due to H<->D transfers.
    //       Should modify precons so size=1 is non-blocking
    level->coarsen(o_r, o_rC);
    parAlmond.Operator(o_rC, o_zC);

    //Add contributions from all patches together
    mesh.ringHalo->Combine(o_zPatch, mesh.Np, ogs_dfloat);

    // Weight by overlap degree, Mr = patchWeight*zPatch
    elliptic.linAlg.amxpy(Ntotal, 1.0, o_patchWeight, o_zPatch, 0.0, o_Mr);

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
  parAlmond(elliptic.platform, settings) {

  //build the one ring mesh
  if (mesh.N>1) {
    meshPatch = mesh.SetupRingPatch();
    ellipticPatch = elliptic.SetupRingPatch(*meshPatch);
    preconPatch = new MultiGridPrecon(*ellipticPatch);
  }

  //build the coarse precon
  int Nf = mesh.N;
  int Nc = 1;  //hard code
  int NpFine   = mesh.Np;
  occa::memory o_weightF = elliptic.o_weight;

  //build mesh and elliptic objects for this degree
  mesh_t &meshC = mesh.SetupNewDegree(Nc);
  elliptic_t &ellipticC = elliptic.SetupNewDegree(meshC);

  //build full A matrix and pass to parAlmond
  parAlmond::parCOO A(elliptic.platform);
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    ellipticC.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    ellipticC.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = elliptic.platform.rank;
  int size = elliptic.platform.size;
  hlong TotalRows = A.globalStarts[size];
  dlong numLocalRows = (dlong) (A.globalStarts[rank+1]-A.globalStarts[rank]);
  dfloat *null = (dfloat *) malloc(numLocalRows*sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1.0/sqrt(TotalRows);

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, elliptic.allNeumann, null,elliptic.allNeumannPenalty);


  // if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
  //   if (parAlmondHandle->numLevels > numMGLevels+1) {
  //     //tell parAlmond to gather when going to the next level
  //     parAlmond::agmgLevel *nextLevel
  //           = (parAlmond::agmgLevel*)parAlmondHandle->levels[numMGLevels+1];

  //     nextLevel->gatherLevel = true;
  //     nextLevel->ogs = ellipticC.ogsMasked;
  //     nextLevel->o_gatherWeight = ellipticC.o_weightG;
  //     nextLevel->Gx = (dfloat*) calloc(nextLevel->R->Ncols,sizeof(dfloat));
  //     nextLevel->Sx = (dfloat*) calloc(meshC.Np*meshC.Nelements,sizeof(dfloat));
  //     nextLevel->o_Gx = elliptic.platform.malloc(nextLevel->R->Ncols*sizeof(dfloat),nextLevel->Gx);
  //     nextLevel->o_Sx = elliptic.platform.malloc(meshC.Np*meshC.Nelements*sizeof(dfloat),nextLevel->Sx);
  //   } else {
  //     //this level is the base
  //     parAlmond::coarseSolver *coarseLevel = parAlmondHandle->coarseLevel;

  //     coarseLevel->gatherLevel = true;
  //     coarseLevel->ogs = ellipticC.ogsMasked;
  //     coarseLevel->Gx = (dfloat*) calloc(coarseLevel->ogs->Ngather,sizeof(dfloat));
  //     coarseLevel->o_Gx = elliptic.platform.malloc(coarseLevel->ogs->Ngather*sizeof(dfloat),coarseLevel->Gx);
  //   }
  // }

  if (mesh.N>1) {
    level = new MGLevel(ellipticC, Nf, NpFine, o_weightF);

    rPatch = (dfloat*) calloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),sizeof(dfloat));
    zPatch = (dfloat*) calloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),sizeof(dfloat));

    dlong Ncols = parAlmond.getNumCols(0);

    rC = (dfloat*) calloc(Ncols,sizeof(dfloat));
    zC = (dfloat*) calloc(Ncols,sizeof(dfloat));

    o_rPatch = elliptic.platform.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), rPatch);
    o_zPatch = elliptic.platform.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), zPatch);

    o_rC = elliptic.platform.malloc(Ncols*sizeof(dfloat), rC);
    o_zC = elliptic.platform.malloc(Ncols*sizeof(dfloat), zC);

    //compute patch overlap weighting
    patchWeight = (dfloat*) malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat));
    for (int i=0;i<mesh.Np*(mesh.Nelements+mesh.totalRingElements);i++)
      patchWeight[i] = 1.0;

    mesh.ringHalo->Combine(patchWeight, mesh.Np, ogs_dfloat);

    //invert
    for (int i=0;i<mesh.Np*(mesh.Nelements+mesh.totalRingElements);i++)
      patchWeight[i] = 1.0/patchWeight[i];

    o_patchWeight = elliptic.platform.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), patchWeight);
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