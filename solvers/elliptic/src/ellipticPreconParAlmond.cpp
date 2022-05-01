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

//AMG preconditioner via parAlmond
void ParAlmondPrecon::Operator(deviceMemory<dfloat>& o_r, deviceMemory<dfloat>& o_Mr) {

  //hand off to parAlmond
  parAlmond.Operator(o_r, o_Mr);

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

ParAlmondPrecon::ParAlmondPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings, elliptic.mesh.comm) {

  //build full A matrix and pass to parAlmond
  if (Comm::World().rank()==0){
    printf("-----------------------------Multigrid AMG Setup--------------------------------------------\n");
  }
  parAlmond::parCOO A(elliptic.platform, elliptic.mesh.comm);
  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    elliptic.BuildOperatorMatrixIpdg(A);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    elliptic.BuildOperatorMatrixContinuous(A);
  }

  //populate null space unit vector
  int rank = elliptic.mesh.rank;
  int size = elliptic.mesh.size;
  hlong TotalRows = A.globalRowStarts[size];
  dlong numLocalRows = static_cast<dlong>(A.globalRowStarts[rank+1]-A.globalRowStarts[rank]);
  memory<dfloat> null(numLocalRows);
  for (dlong i=0;i<numLocalRows;i++) {
    null[i] = 1.0/sqrt(TotalRows);
  }

  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);

  parAlmond.Report();

  //The csr matrix at the top level of parAlmond may have a larger
  // halo region than the matrix free kernel. Adjust if necessary
  dlong parAlmondNrows = parAlmond.getNumRows(0);
  dlong parAlmondNcols = parAlmond.getNumCols(0);
  dlong parAlmondNhalo = parAlmondNcols - parAlmondNrows;
  _elliptic.Nhalo = std::max(_elliptic.Nhalo, parAlmondNhalo);
}
