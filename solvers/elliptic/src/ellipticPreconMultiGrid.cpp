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


// Matrix-free p-Multigrid levels followed by AMG
void MultiGridPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  //just pass to parAlmond
  parAlmond::Precon(parAlmondHandle, o_Mr, o_r);

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

MultiGridPrecon::MultiGridPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings) {

  //initialize parAlmond
  parAlmondHandle = parAlmond::Init(mesh.device, mesh.comm, settings);
  parAlmond::multigridLevel **levels = parAlmondHandle->levels;

  int Nf = mesh.N;
  int Nc = Nf;
  int NpFine   = mesh.Np;
  int NpCoarse = mesh.Np;
  occa::memory o_weightF = elliptic.o_weight;

  NpMGlevels=0;

  while(true) {
    //build mesh and elliptic objects for this degree
    mesh_t &meshC = mesh.SetupNewDegree(Nc);
    elliptic_t &ellipticC = elliptic.SetupNewDegree(meshC);
    NpMGlevels++;

    if (Nc==1) { //base p-MG level
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
      break;

    } else {
      //make a multigrid level
      int numMGLevels = parAlmondHandle->numLevels;
      levels[numMGLevels] = new MGLevel(ellipticC, numMGLevels, Nf, NpFine, o_weightF,
                                        parAlmondHandle->ktype, parAlmondHandle->ctype);
      parAlmondHandle->numLevels++;
    }

    //find the degree of the next level
    Nf = Nc;
    NpFine = meshC.Np;
    o_weightF = ellipticC.o_weight; //save previous weights
    if (settings.compareSetting("MULTIGRID COARSENING","ALLDEGREES")) {
      Nc = Nf-1;
    } else if (settings.compareSetting("MULTIGRID COARSENING","HALFDEGREES")) {
      Nc = (Nf+1)/2;
    } else { //default "HALFDOFS"
      // pick the degrees so the dofs of each level halfs (roughly)
      while (NpCoarse > NpFine/2) {
        Nc--;
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
      }
    }
  }

  //report
  if (settings.compareSetting("VERBOSE","TRUE")) {
    //This setup can be called by many subcommunicators, so only
    // print on the global root.
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank==0) { //report the upper multigrid levels
      printf("------------------Multigrid Report----------------------------------------\n");
      printf("--------------------------------------------------------------------------\n");
      printf("level|    Type    |    dimension   |   nnz per row   |   Smoother        |\n");
      printf("     |            |  (min,max,avg) |  (min,max,avg)  |                   |\n");
      printf("--------------------------------------------------------------------------\n");
    }

    for(int lev=0; lev<parAlmondHandle->numLevels; lev++) {
      if(rank==0) {printf(" %3d ", lev);fflush(stdout);}
      levels[lev]->Report();
    }

    if (rank==0)
      printf("--------------------------------------------------------------------------\n");
  }
}

MultiGridPrecon::~MultiGridPrecon() {

  for (int i=1;i<NpMGlevels;i++) {
    MGLevel *level = (MGLevel *) parAlmondHandle->levels[i];
    delete &(level->elliptic);
    if (level->mesh.ogs) level->mesh.ogs->Free();
  }

  if (parAlmondHandle) parAlmond::Free(parAlmondHandle);
}