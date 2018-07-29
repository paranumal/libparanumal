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

void boltzmannReportAddaptive2D(bns_t *bns, dfloat t, setupAide &options){

  mesh2D *mesh = bns->mesh; 
  
  // copy data back to host
  bns->o_q.copyTo(bns->q);


  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0){
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
    printf("t: %g ramp: %g drampdt: %g\n", t, ramp, drampdt);
  }
  
  
  //
  #if 1
  boltzmannForces2D(bns,t,options);
  #endif


 
  if(options.compareArgs("ABSORBING LAYER","PML")){ 

    if(options.compareArgs("OUTPUT FILE FORMAT","VTU")){ 
    
      
    char fname[BUFSIZ];
    // sprintf(fname, "/scratch/boltzmannInclined/foo_pml_%.0f_%04d_%04d.vtu", bns->Re, rank, tstep/bns->errorStep);
    sprintf(fname, "foo_pml_%.0f_%04d_%04d.vtu", bns->Re, rank, bns->frame++);
    boltzmannPlotVTU2D(bns, fname);
   }


  if(options.compareArgs("OUTPUT FILE FORMAT","TEC")){ 
    //boltzmannComputeVorticity2D(mesh, mesh->q,5, mesh->Nfields);
    char fname[BUFSIZ];
    sprintf(fname, "foo_v2_%04d.dat",rank);
    // boltzmannPlotTEC2D(mesh, fname, tstep/mesh->errorStep);
    boltzmannPlotTEC2D(bns, fname, t);
  }


  }
  else{

   if(options.compareArgs("OUTPUT FILE FORMAT","VTU")){ 
    //boltzmannCouetteError2D(mesh, t);
    // compute vorticity
    //boltzmannComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    // int fld = 1;
    char fname[BUFSIZ];
    sprintf(fname, "foo2_%04d_%04d.vtu",rank, bns->frame++);
    boltzmannPlotVTU2D(bns, fname);
  }

  
  if(options.compareArgs("OUTPUT FILE FORMAT","TEC")){ 
    char fname[BUFSIZ];
    sprintf(fname, "foo_%04d.vtu",rank);
    boltzmannPlotTEC2D(bns, fname, t);
  }    
  }
  
}
