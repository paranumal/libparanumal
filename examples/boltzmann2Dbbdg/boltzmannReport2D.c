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

#include "boltzmann2D.h"

void boltzmannReport2D(mesh2D *mesh, int tstep){

  dfloat t = (tstep+1)*mesh->dt;
  
  // copy data back to host
  mesh->o_q.copyTo(mesh->q);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
    printf("t: %g ramp: %g drampdt: %g\n", t, ramp, drampdt);
  }
  
  // do error stuff on host
  boltzmannError2D(mesh, t);
  
  // compute vorticity
  boltzmannComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
  
  // output field files
  int fld = 0;
  char fname[BUFSIZ];
  sprintf(fname, "foo_T%04d", tstep/mesh->errorStep);
  meshPlotVTU2D(mesh, fname, fld);
}
