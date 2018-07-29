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

#include "insQuad2D.h"

void insReportQuad2D(ins_t *ins, int tstep, char *options){
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh2D *mesh = ins->mesh;

  dfloat t = (tstep)*ins->dt;
  
  dlong offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  ins->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_D,
                       offset,
                       ins->o_U,
                       ins->o_V,
                       ins->o_Vort);

  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_D,
                             offset,
                             ins->o_U,
                             ins->o_V,
                             ins->o_Div);

  // gather-scatter
  ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, ins->o_Vort, dfloatString, "add");  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vort, ins->o_Vort);

  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V); 
  ins->o_P.copyTo(ins->P);

  ins->o_Vort.copyTo(ins->Vort);
  ins->o_Div.copyTo(ins->Div);
  
  // do error stuff on host
  insErrorQuad2D(ins, t, options);
 
  if(strstr(options, "VTU")){ 
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins2D/foo_%04d_%04d.vtu",rank, tstep/ins->errorStep);
    sprintf(fname, "foo_%04d_%04d.vtu",rank, tstep/ins->errorStep);

    insPlotVTUQuad2D(ins, fname);
  } 
}

