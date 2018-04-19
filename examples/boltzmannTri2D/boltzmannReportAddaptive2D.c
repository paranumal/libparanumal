#include "boltzmann2D.h"

void boltzmannReportAddaptive2D(bns_t *bns, dfloat t, char *options){

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


  // do error stuff on host
  boltzmannError2D(bns, t, options);




  if(strstr(options, "PML")){ 

    if(strstr(options, "VTU")){ 
    
    #if 0
    char fname[BUFSIZ];
    sprintf(fname, "Allfields_Obl_mpml__t3_%04d_%04d.vtu", rank, tstep/mesh->errorStep);
    boltzmannPlotVTUField2D(mesh, fname);
    #endif
      
    char fname[BUFSIZ];
    // sprintf(fname, "/scratch/boltzmannInclined/foo_pml_%.0f_%04d_%04d.vtu", bns->Re, rank, tstep/bns->errorStep);
    sprintf(fname, "foo_pml_%.0f_%04d_%04d.vtu", bns->Re, rank, bns->frame++);
    boltzmannPlotVTU2D(bns, fname);
   }


  if(strstr(options, "TEC")){ 
    //boltzmannComputeVorticity2D(mesh, mesh->q,5, mesh->Nfields);
    char fname[BUFSIZ];
    sprintf(fname, "foo_v2_%04d.dat",rank);
    // boltzmannPlotTEC2D(mesh, fname, tstep/mesh->errorStep);
    boltzmannPlotTEC2D(bns, fname, t);
  }


  }
  else{

   if(strstr(options, "VTU")){ 
    //boltzmannCouetteError2D(mesh, t);
    // compute vorticity
    //boltzmannComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    // int fld = 1;
    char fname[BUFSIZ];
    sprintf(fname, "foo2_%04d_%04d.vtu",rank, bns->frame++);
    boltzmannPlotVTU2D(bns, fname);
  }

  
  if(strstr(options, "TEC")){ 
    char fname[BUFSIZ];
    sprintf(fname, "foo_%04d.vtu",rank);
    boltzmannPlotTEC2D(bns, fname, bns->frame++);
  }    
  }
  
}
