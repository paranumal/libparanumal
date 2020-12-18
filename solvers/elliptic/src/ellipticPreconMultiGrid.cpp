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

  if (gather==true) { //gather before passing to parAlmond
    elliptic.ogsMasked->Gather(o_rhsG, o_r, ogs_dfloat, ogs_add, ogs_notrans);
    parAlmond.Operator(o_rhsG, o_xG);
    elliptic.ogsMasked->Scatter(o_Mr, o_xG, ogs_dfloat, ogs_add, ogs_notrans);
  } else {
    //just pass to parAlmond
    parAlmond.Operator(o_r, o_Mr);

    if (elliptic.disc_c0) {
      dlong Ntotal = elliptic.mesh.Nelements*elliptic.mesh.Np;
      elliptic.ogsMasked->GatherScatter(o_Mr, ogs_dfloat, ogs_add, ogs_sym);
      elliptic.linAlg.amx(Ntotal, (dfloat)1.0, elliptic.o_weight, o_Mr);
    }
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

MultiGridPrecon::MultiGridPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings, mesh.comm) {

  int Nf = mesh.N;
  int Nc = Nf;
  int NpFine   = mesh.Np;
  int NpCoarse = mesh.Np;

  MGLevel* prevLevel=nullptr;
  MGLevel* currLevel=nullptr;

  int minN = 1; // degree of coarsest mesh
  
  printf("allNeumannPenalty=%g\n", elliptic.allNeumannPenalty);
  
  while(Nc>minN) {
    //build mesh and elliptic objects for this degree
    mesh_t &meshF = mesh.SetupNewDegree(Nf);
    elliptic_t &ellipticF = elliptic.SetupNewDegree(meshF);

    //share masking data with previous MG level
    if (prevLevel) {
      prevLevel->Nmasked = ellipticF.Nmasked;
      prevLevel->o_maskIds = ellipticF.o_maskIds;
      prevLevel->ogsMasked = ellipticF.ogsMasked;
    }

    //find the degree of the next level
    if (settings.compareSetting("MULTIGRID COARSENING","ALLDEGREES")) {
      Nc = mymax(1,Nf-1);
    } else if (settings.compareSetting("MULTIGRID COARSENING","HALFDEGREES")) {
      Nc = mymax(1,(Nf+1)/2);
    } else { //default "HALFDOFS"
      // pick the degrees so the dofs of each level halfs (roughly)
      while (NpCoarse > NpFine/2 && Nc>1) {
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
    
    //make a multigrid level
    currLevel = new MGLevel(ellipticF, Nc, NpCoarse);
    parAlmond.AddLevel(currLevel);

    Nf = Nc;
    NpFine = NpCoarse;
    prevLevel = currLevel;
  }

  //build matrix at degree 1
  mesh_t &meshF = mesh.SetupNewDegree(minN);
  elliptic_t &ellipticF = elliptic.SetupNewDegree(meshF);

  //build full A matrix and pass to parAlmond
  parAlmond::parCOO A(elliptic.platform, mesh.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    ellipticF.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    ellipticF.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = mesh.rank;
  int size = mesh.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = (dlong) (A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);
  dfloat *null = (dfloat *) malloc(numLocalRows*sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1.0/sqrt(TotalRows);

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);
  free(null);

  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    if (mesh.N>minN) {
      //tell the last pMG level to gather after coarsening
      prevLevel->gatherLevel = true;
      prevLevel->ogsMasked = ellipticF.ogsMasked;

      dfloat *dummy = (dfloat *) calloc(meshF.Np*meshF.Nelements,sizeof(dfloat));
      prevLevel->o_SX = elliptic.platform.malloc(meshF.Np*meshF.Nelements*sizeof(dfloat), dummy);
      prevLevel->o_GX = elliptic.platform.malloc(meshF.Np*meshF.Nelements*sizeof(dfloat), dummy);
      free(dummy);
    } else {
      //gather before passing to parAlmond
      gather=true;
      dlong Ncols = parAlmond.getNumCols(0);
      o_rhsG = elliptic.platform.malloc(Ncols*sizeof(dfloat));
      o_xG   = elliptic.platform.malloc(Ncols*sizeof(dfloat));
    }
  }

  //report
  parAlmond.Report();
}
