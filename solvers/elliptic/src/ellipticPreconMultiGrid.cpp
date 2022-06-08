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


// Matrix-free p-Multigrid levels followed by AMG
void MultiGridPrecon::Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr) {

  //just pass to parAlmond
  parAlmond.Operator(o_r, o_Mr);

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

  while(Nc>1) {
    if (Comm::World().rank()==0){
      printf("-----------------------------Multigrid pMG Degree %2d----------------------------------------\n", Nc);
    }
    //build mesh and elliptic objects for this degree
    mesh_t meshF = mesh.SetupNewDegree(Nf);
    elliptic_t ellipticF = elliptic.SetupNewDegree(meshF);

    //share masking data with previous MG level
    if (parAlmond.NumLevels()>0) {
      MGLevel& prevLevel = parAlmond.GetLevel<MGLevel>(parAlmond.NumLevels()-1);
      prevLevel.meshC = meshF;
      prevLevel.ellipticC = ellipticF;
    }

    //find the degree of the next level
    if (settings.compareSetting("MULTIGRID COARSENING","ALLDEGREES")) {
      Nc = Nf-1;
    } else if (settings.compareSetting("MULTIGRID COARSENING","HALFDEGREES")) {
      Nc = std::max(1,(Nf+1)/2);
    } else { //default "HALFDOFS"
      // pick the degrees so the dofs of each level halfs (roughly)
      while (NpCoarse > NpFine/2 && Nc>1) {
        Nc--;
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
      }
    }

    //set Npcoarse
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

    dlong Nrows, Ncols;
    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      Nrows = ellipticF.ogsMasked.Ngather;
      Ncols = Nrows + ellipticF.gHalo.Nhalo;
    } else {
      Nrows = meshF.Nelements*meshF.Np;
      Ncols = Nrows + meshF.totalHaloPairs*mesh.Np;
    }

    //Add a multigrid level
    parAlmond.AddLevel<MGLevel>(ellipticF, Nrows, Ncols, Nc, NpCoarse);

    Nf = Nc;
    NpFine = NpCoarse;
  }

  //build matrix at degree 1
  if (Comm::World().rank()==0){
    printf("-----------------------------Multigrid pMG Degree  1----------------------------------------\n");
  }
  mesh_t meshF = mesh.SetupNewDegree(1);
  elliptic_t ellipticF = elliptic.SetupNewDegree(meshF);

  //share masking data with previous MG level
  if (parAlmond.NumLevels()>0) {
    MGLevel& prevLevel = parAlmond.GetLevel<MGLevel>(parAlmond.NumLevels()-1);
    prevLevel.meshC = meshF;
    prevLevel.ellipticC = ellipticF;
  }

  //build full A matrix and pass to parAlmond
  if (Comm::World().rank()==0){
    printf("-----------------------------Multigrid AMG Setup--------------------------------------------\n");
  }
  parAlmond::parCOO A(elliptic.platform, mesh.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG"))
    ellipticF.BuildOperatorMatrixIpdg(A);
  else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS"))
    ellipticF.BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = mesh.rank;
  int size = mesh.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = static_cast<dlong>(A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);

  memory<dfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  //set up AMG levels (treating the N=1 level as a matrix level)
  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);

  //report
  parAlmond.Report();
}
