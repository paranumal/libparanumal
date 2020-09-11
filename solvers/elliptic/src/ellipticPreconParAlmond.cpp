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

//AMG preconditioner via parAlmond
void ParAlmondPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  if (elliptic.disc_c0) {
    elliptic.ogsMasked->Gather(o_rhsG, o_r, ogs_dfloat, ogs_add, ogs_notrans);
    parAlmond.Operator(o_rhsG, o_xG);
    elliptic.ogsMasked->Scatter(o_Mr, o_xG, ogs_dfloat, ogs_add, ogs_notrans);
  } else {
    //IPDG
    parAlmond.Operator(o_r, o_Mr);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

ParAlmondPrecon::ParAlmondPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings, elliptic.mesh.comm) {

  //build full A matrix and pass to parAlmond
  parAlmond::parCOO A(elliptic.platform, elliptic.mesh.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    elliptic.BuildOperatorMatrixIpdg(A);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    elliptic.BuildOperatorMatrixContinuous(A);
  }

  //populate null space unit vector
  int rank = elliptic.mesh.rank;
  int size = elliptic.mesh.size;
  hlong TotalRows = A.globalStarts[size];
  dlong numLocalRows = (dlong) (A.globalStarts[rank+1]-A.globalStarts[rank]);
  dfloat *null = (dfloat *) malloc(numLocalRows*sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1.0/sqrt(TotalRows);

  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);
  free(null);

  parAlmond.Report();

  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    //make buffers to gather this level before passing to parAlmond
    dlong Ncols = parAlmond.getNumCols(0);

    rhsG = (dfloat*) calloc(Ncols,sizeof(dfloat));
    xG   = (dfloat*) calloc(Ncols,sizeof(dfloat));
    o_rhsG = elliptic.platform.malloc(Ncols*sizeof(dfloat));
    o_xG   = elliptic.platform.malloc(Ncols*sizeof(dfloat));
  }
}

ParAlmondPrecon::~ParAlmondPrecon() {
  if (rhsG) free(rhsG);
  if (xG) free(xG);
}