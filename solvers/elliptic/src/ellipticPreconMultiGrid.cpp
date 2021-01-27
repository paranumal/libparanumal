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

void occaConvertType(platform_t &platform, dlong N, occa::memory &o_x, string &gfloatString){


  if(gfloatString==dfloatString) return;

  // bring o_x data back to host
  dfloat *h_x = (dfloat*) calloc(N, sizeof(dfloat));
  o_x.copyTo(h_x);

  // free device storage
  o_x.free();

  // type convert on host and copy back to device
  if(gfloatString=="float"){
    float *x = (float*) calloc(N, sizeof(float));
    for(int n=0;n<N;++n){
      x[n] = (float)h_x[n];
    }
    o_x = platform.device.malloc(N*sizeof(float), x);
    free(x);
  }

  if(gfloatString=="double"){
    double *x = (double*) calloc(N, sizeof(double));
    for(int n=0;n<N;++n){
      x[n] = (double)h_x[n];
    }
    o_x = platform.device.malloc(N*sizeof(double), x);
    free(x);
  }

  free(h_x);
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

  while(Nc>1) {
    //build mesh and elliptic objects for this degree
    mesh_t &meshF = mesh.SetupNewDegree(Nf);

    // TW: rewrite o_ggeo here ? and set a variable o_ggeoType in mesh ?
    if(Nf<mesh.N){ //
      meshF.gfloatString = dfloatString;
      settings.getSetting("MULTIGRID GEOFAC TYPE", meshF.gfloatString);

      dlong M = 0;
      if(meshF.elementType==TRIANGLES || meshF.elementType==TETRAHEDRA)
	M = meshF.Nelements*meshF.Nggeo;
      else
	M = meshF.Nelements*meshF.Np*meshF.Nggeo;

      occaConvertType(meshF.platform, M, meshF.o_ggeo, meshF.gfloatString);
    }
    
    elliptic_t &ellipticF = elliptic.SetupNewDegree(meshF);

    //share masking data with previous MG level
    if (prevLevel) {
      prevLevel->meshC = &meshF;
      prevLevel->ogsMaskedC = ellipticF.ogsMasked;
    }

    //find the degree of the next level
    if (settings.compareSetting("MULTIGRID COARSENING","ALLDEGREES")) {
      Nc = Nf-1;
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

    //set Npcoarse
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

    dlong Nrows, Ncols;
    if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
      Nrows = ellipticF.ogsMasked->Ngather;
      Ncols = Nrows + ellipticF.ogsMasked->NgatherHalo;
    } else {
      Nrows = meshF.Nelements*meshF.Np;
      Ncols = Nrows + meshF.totalHaloPairs*mesh.Np;
    }

    //make a multigrid level
    currLevel = new MGLevel(ellipticF, Nrows, Ncols, Nc, NpCoarse);
    parAlmond.AddLevel(currLevel);

    Nf = Nc;
    NpFine = NpCoarse;
    prevLevel = currLevel;
  }

  //build matrix at degree 1
  mesh_t &meshF = mesh.SetupNewDegree(1);
  elliptic_t &ellipticF = elliptic.SetupNewDegree(meshF);

  //share masking data with previous MG level
  if (prevLevel) {
    prevLevel->meshC = &meshF;
    prevLevel->ogsMaskedC = ellipticF.ogsMasked;
  }

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

  //report
  parAlmond.Report();
}
