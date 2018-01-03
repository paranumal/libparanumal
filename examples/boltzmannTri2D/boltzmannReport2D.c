#include "boltzmann2D.h"

void boltzmannReport2D(mesh2D *mesh, iint tstep, char *options){
  
  dfloat t = 0.f; 

  if(strstr(options,"MRAB") || strstr(options, "MRSAAB"))
   t = mesh->dt*(tstep+1)*pow(2,(mesh->MRABNlevels-1));
  else
   t = (tstep+1)*mesh->dt;
  
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
  
  

  if(strstr(options, "PML")){ 
    // do error stuff on host
    boltzmannError2D(mesh, t, options);

   
    if(strstr(options, "VTU")){ 
    // compute vorticity
    //boltzmannComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    iint fld = 1;
    char fname[BUFSIZ];
    sprintf(fname, "foo_%04d_%04d.vtu", rank, tstep/mesh->errorStep);
    boltzmannPlotVTU2D(mesh, fname);
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
    sprintf(fname, "foo_%04d_%04d.vtu",rank, tstep/mesh->errorStep);
    boltzmannPlotVTU2D(mesh, fname);
  }
  
   // do error stuff on host
    boltzmannError2D(mesh, t, options);
    
  }
  
}
