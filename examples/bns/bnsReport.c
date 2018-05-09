#include "bns.h"

void bnsReport(bns_t *bns,  int tstep, setupAide &options){

  mesh_t *mesh = bns->mesh; 
  
  dfloat t = 0.f; 

  if(options.compareArgs("TIME INTEGRATOR","MRSAAB"))
   t = bns->startTime + bns->dt*tstep*pow(2,(mesh->MRABNlevels-1));
  else
   t = bns->startTime + tstep*bns->dt;
  
  // copy data back to host
  bns->o_q.copyTo(bns->q);


  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0){
    dfloat ramp, drampdt;
    bnsRampFunction(t, &ramp, &drampdt);
    printf("t: %g ramp: %g drampdt: %g\n", t, ramp, drampdt);
  }
  
  
  // not implemented yet
  #if 0 
    bnsForces(bns,t,options);
  #endif


  if(options.compareArgs("OUTPUT FILE FORMAT","VTU")){ 
    char fname[BUFSIZ];
    // sprintf(fname, "/scratch/boltzmannInclined/foo_pml_%.0f_%04d_%04d.vtu", bns->Re, rank, tstep/bns->errorStep);
    sprintf(fname, "foo_pml_%04d_%04d.vtu",rank, tstep/bns->errorStep);
    bnsPlotVTU(bns, fname);
  }


  if(options.compareArgs("OUTPUT FILE FORMAT","TEC")){ 
    // //boltzmannComputeVorticity2D(mesh, mesh->q,5, mesh->Nfields);
    // char fname[BUFSIZ];
    // sprintf(fname, "foo_v2_%04d.dat",rank);
    // bnsPlotTEC(bns, fname, t);
  }
  
  
}
