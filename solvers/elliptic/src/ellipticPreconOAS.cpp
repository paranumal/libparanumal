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
    parAlmond::multigridLevel *level = parAlmondHandle->levels[0];
    level->coarsen(o_r, o_rC);
    parAlmond::Precon(parAlmondHandle, o_zC, o_rC);

    //Add contributions from all patches together
    mesh.ringHalo->Combine(o_zPatch, mesh.Np, ogs_dfloat);

    // Weight by overlap degree, Mr = patchWeight*zPatch
    elliptic.linAlg.amxpy(Ntotal, 1.0, o_patchWeight, o_zPatch, 0.0, o_Mr);

    // Add prologatated coarse solution
    level->prolongate(o_zC, o_Mr);
  } else {
    //if N=1 just call the coarse solver
    parAlmond::Precon(parAlmondHandle, o_Mr, o_r);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

OASPrecon::OASPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings) {

  //initialize parAlmond
  parAlmondHandle = parAlmond::Init(mesh.device, mesh.comm, settings);
  parAlmond::multigridLevel **levels = parAlmondHandle->levels;

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
  parAlmond::parCOO A;
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    ellipticC.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    ellipticC.BuildOperatorMatrixContinuous(A);

  int numMGLevels = parAlmondHandle->numLevels;

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond::AMGSetup(parAlmondHandle, A,
                  elliptic.allNeumann, elliptic.allNeumannPenalty);

  //overwrite the finest AMG level with the degree 1 matrix free level
  delete levels[numMGLevels];

  levels[numMGLevels] = new MGLevel(ellipticC, numMGLevels, Nf, NpFine, o_weightF,
                                    parAlmondHandle->ktype, parAlmondHandle->ctype);

  if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    if (parAlmondHandle->numLevels > numMGLevels+1) {
      //tell parAlmond to gather when going to the next level
      parAlmond::agmgLevel *nextLevel
            = (parAlmond::agmgLevel*)parAlmondHandle->levels[numMGLevels+1];

      nextLevel->gatherLevel = true;
      nextLevel->ogs = ellipticC.ogsMasked;
      nextLevel->o_gatherWeight = ellipticC.o_weightG;
      nextLevel->Gx = (dfloat*) calloc(nextLevel->R->Ncols,sizeof(dfloat));
      nextLevel->Sx = (dfloat*) calloc(meshC.Np*meshC.Nelements,sizeof(dfloat));
      nextLevel->o_Gx = meshC.device.malloc(nextLevel->R->Ncols*sizeof(dfloat),nextLevel->Gx);
      nextLevel->o_Sx = meshC.device.malloc(meshC.Np*meshC.Nelements*sizeof(dfloat),nextLevel->Sx);
    } else {
      //this level is the base
      parAlmond::coarseSolver *coarseLevel = parAlmondHandle->coarseLevel;

      coarseLevel->gatherLevel = true;
      coarseLevel->ogs = ellipticC.ogsMasked;
      coarseLevel->Gx = (dfloat*) calloc(coarseLevel->ogs->Ngather,sizeof(dfloat));
      coarseLevel->o_Gx = meshC.device.malloc(coarseLevel->ogs->Ngather*sizeof(dfloat),coarseLevel->Gx);
    }
  }

  if (mesh.N>1) {
    rPatch = (dfloat*) calloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),sizeof(dfloat));
    zPatch = (dfloat*) calloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements),sizeof(dfloat));

    rC = (dfloat*) calloc(levels[0]->Ncols,sizeof(dfloat));
    zC = (dfloat*) calloc(levels[0]->Ncols,sizeof(dfloat));

    o_rPatch = mesh.device.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), rPatch);
    o_zPatch = mesh.device.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), zPatch);

    o_rC = mesh.device.malloc(levels[0]->Ncols*sizeof(dfloat), rC);
    o_zC = mesh.device.malloc(levels[0]->Ncols*sizeof(dfloat), zC);

    //compute patch overlap weighting
    patchWeight = (dfloat*) malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat));
    for (int i=0;i<mesh.Np*(mesh.Nelements+mesh.totalRingElements);i++)
      patchWeight[i] = 1.0;

    mesh.ringHalo->Combine(patchWeight, mesh.Np, ogs_dfloat);

    //invert
    for (int i=0;i<mesh.Np*(mesh.Nelements+mesh.totalRingElements);i++)
      patchWeight[i] = 1.0/patchWeight[i];

    o_patchWeight = mesh.device.malloc(mesh.Np*(mesh.Nelements+mesh.totalRingElements)*sizeof(dfloat), patchWeight);
  }

  //report
  if (settings.compareSetting("VERBOSE","TRUE")) {
    if (mesh.rank==0) { //report the upper multigrid levels
      printf("------------------Multigrid Report----------------------------------------\n");
      printf("--------------------------------------------------------------------------\n");
      printf("level|    Type    |    dimension   |   nnz per row   |   Smoother        |\n");
      printf("     |            |  (min,max,avg) |  (min,max,avg)  |                   |\n");
      printf("--------------------------------------------------------------------------\n");
    }

    for(int lev=0; lev<parAlmondHandle->numLevels; lev++) {
      if(mesh.rank==0) {printf(" %3d ", lev);fflush(stdout);}
      levels[lev]->Report();
    }

    if (mesh.rank==0)
      printf("--------------------------------------------------------------------------\n");
  }
}

OASPrecon::~OASPrecon() {
  if (mesh.N>1) {
    delete preconPatch;
    if (mesh.size>1) delete ellipticPatch;
    if (mesh.size>1) delete meshPatch;

    MGLevel *level = (MGLevel *) parAlmondHandle->levels[0];
    delete &(level->elliptic);
    if (level->mesh.ogs) level->mesh.ogs->Free();
  }
  if (parAlmondHandle) parAlmond::Free(parAlmondHandle);
}