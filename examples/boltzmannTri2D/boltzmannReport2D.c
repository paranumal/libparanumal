#include "boltzmann2D.h"

void boltzmannReport2D(mesh2D *mesh, iint tstep, char *options){
  
  dfloat t = 0.f; 

  if(strstr(options,"MRAB") || strstr(options, "MRSAAB"))
   t = mesh->startTime + mesh->dt*tstep*pow(2,(mesh->MRABNlevels-1));
  else
   t = mesh->startTime+ tstep*mesh->dt;
  
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
  
  
  //
  #if 1
  boltzmannForces2D(mesh,t,options);
  #endif


  // do error stuff on host
  boltzmannError2D(mesh, t, options);




  if(strstr(options, "PML")){ 

    if(strstr(options, "VTU")){ 
    
    #if 0
    char fname[BUFSIZ];
    sprintf(fname, "Allfields_Obl_mpml__t3_%04d_%04d.vtu", rank, tstep/mesh->errorStep);
    boltzmannPlotVTUField2D(mesh, fname);
    #endif
      
    char fname[BUFSIZ];
    sprintf(fname, "foo_pml_%.0f_%04d_%04d.vtu", mesh->Re, rank, tstep/mesh->errorStep);
    boltzmannPlotVTU2D(mesh, fname);
   }


  if(strstr(options, "TEC")){ 
    //boltzmannComputeVorticity2D(mesh, mesh->q,5, mesh->Nfields);
    char fname[BUFSIZ];
    sprintf(fname, "foo_v2_%04d.dat",rank);
    // boltzmannPlotTEC2D(mesh, fname, tstep/mesh->errorStep);
    boltzmannPlotTEC2D(mesh, fname, t);
  }


  }
  else{

   if(strstr(options, "VTU")){ 
    //boltzmannCouetteError2D(mesh, t);
    // compute vorticity
    //boltzmannComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    // iint fld = 1;
    char fname[BUFSIZ];
    sprintf(fname, "foo2_%04d_%04d.vtu",rank, tstep/mesh->errorStep);
    boltzmannPlotVTU2D(mesh, fname);
  }

  
  if(strstr(options, "TEC")){ 
    char fname[BUFSIZ];
    sprintf(fname, "foo_%04d.vtu",rank);
    boltzmannPlotTEC2D(mesh, fname, tstep/mesh->errorStep);
  }    
  }
  
}
