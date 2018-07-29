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

#include "mns.h"

void mnsReport(mns_t *mns, dfloat time, int tstep){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = mns->mesh;

  // mns->vorticityKernel(mesh->Nelements,
  //                      mesh->o_vgeo,
  //                      mesh->o_Dmatrices,
  //                      mns->fieldOffset,
  //                      mns->o_U,
  //                      mns->o_Vort);

  // mns->divergenceVolumeKernel(mesh->Nelements,
  //                            mesh->o_vgeo,
  //                            mesh->o_Dmatrices,
  //                            mns->fieldOffset,
  //                            mns->o_U,
  //                            mns->o_Div);

  // gather-scatter
  // ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_Vort, dfloatString, "add");  
  // ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vort, ins->o_Vort);
  
  // copy data back to host
  mns->o_U.copyTo(mns->U);
  mns->o_P.copyTo(mns->P);
  mns->o_Phi.copyTo(mns->Phi);
  // mns->o_Vort.copyTo(mns->Vort);
  // mns->o_Div.copyTo(mns->Div);
  
  mns->o_GPhi.copyTo(mns->GPhi);
  mns->o_SPhi.copyTo(mns->SPhi);
  // do error stuff on host
  mnsError(mns, time);

  if(mns->options.compareArgs("OUTPUT TYPE","VTU")){ 
    // output field files
    char fname[BUFSIZ];
    string outName;
    mns->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(),rank, tstep/mns->outputStep);

    mnsPlotVTU(mns, fname);
  } 
}

