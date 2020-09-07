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
  parAlmond.Operator(o_r, o_Mr);

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

MultiGridPrecon::MultiGridPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings) {

  int Nf = mesh.N;
  int Nc = Nf;
  int NpFine   = mesh.Np;
  int NpCoarse = mesh.Np;
  occa::memory o_weightF = elliptic.o_weight;

  while(true) {
    //build mesh and elliptic objects for this degree
    mesh_t &meshC = mesh.SetupNewDegree(Nc);
    elliptic_t &ellipticC = elliptic.SetupNewDegree(meshC);

    if (Nc==1) { //base p-MG level
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
      parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);
      free(null);

      // int numMGLevels = parAlmondHandle->numLevels;
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
      break;

    } else {
      //make a multigrid level
      MGLevel* level = new MGLevel(ellipticC, Nf, NpFine, o_weightF);
      parAlmond.AddLevel(level);
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
  parAlmond.Report();
}
