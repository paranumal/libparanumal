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
ParAlmondPrecon::ParAlmondPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), settings(_elliptic.settings) {

  mesh_t& mesh = elliptic.mesh;

  //build full A matrix and pass to parAlmond
  parAlmond::parCOO A;
  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {
    elliptic.BuildOperatorMatrixIpdg(A);
  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    elliptic.BuildOperatorMatrixContinuous(A);
  }

  parAlmondHandle = parAlmond::Init(mesh.device, mesh.comm, settings);
  parAlmond::AMGSetup(parAlmondHandle, A,
                      elliptic.allNeumann, elliptic.allNeumannPenalty);

  parAlmond::Report(parAlmondHandle);

  if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    //make buffers to gather this level before passing to parAlmond
    parAlmond::multigridLevel *baseLevel = parAlmondHandle->levels[0];

    rhsG = (dfloat*) calloc(baseLevel->Ncols,sizeof(dfloat));
    xG   = (dfloat*) calloc(baseLevel->Ncols,sizeof(dfloat));
    o_rhsG = mesh.device.malloc(baseLevel->Ncols*sizeof(dfloat));
    o_xG   = mesh.device.malloc(baseLevel->Ncols*sizeof(dfloat));
  }
}

void ParAlmondPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  if (settings.compareSetting("DISCRETIZATION", "IPDG")) {

    parAlmond::Precon(parAlmondHandle, o_Mr, o_r);

  } else if (settings.compareSetting("DISCRETIZATION", "CONTINUOUS")) {
    elliptic.ogsMasked->Gather(o_rhsG, o_r, ogs_dfloat, ogs_add, ogs_notrans);
    parAlmond::Precon(parAlmondHandle, o_xG, o_rhsG);
    elliptic.ogsMasked->Scatter(o_Mr, o_xG, ogs_dfloat, ogs_add, ogs_notrans);
  }

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

ParAlmondPrecon::~ParAlmondPrecon() {
  if (parAlmondHandle) parAlmond::Free(parAlmondHandle);
}