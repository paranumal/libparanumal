#include "boltzmann3D.h"

void boltzmannReport3D(mesh3D *mesh, iint tstep, char *options){

  dfloat t = (tstep+1)*mesh->dt;
  
  // copy data back to host
  mesh->o_q.copyTo(mesh->q);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    dfloat ramp, drampdt;
    boltzmannRampFunction3D(t, &ramp, &drampdt);
    printf("t: %g ramp: %g drampdt: %g\n", t, ramp, drampdt);
  }
  
  

  if(strstr(options, "PML")){ 
    // do error stuff on host
    boltzmannError3D(mesh, t, options);
    if(strstr(options, "VTU")){ 

    //iint fld = 0;  
    // compute vorticity
    //boltzmannComputeVorticity3D(mesh, mesh->q, fld, mesh->Nfields);
    // output field files
    
    char fname[BUFSIZ];
    sprintf(fname, "fooT_%04d", tstep/mesh->errorStep);
    boltzmannPlotVTU3D(mesh, fname);
    }
  }
  else{
    // do error stuff on host
    boltzmannError3D(mesh, t, options);

   if(strstr(options, "VTU")){ 
    //boltzmannCouetteError2D(mesh, t);
    // compute vorticity
    //iint fld = 0;
    //boltzmannComputeVorticity3D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    
    char fname[BUFSIZ];
    sprintf(fname, "fooT_%04d", tstep/mesh->errorStep);
    boltzmannPlotVTU3D(mesh, fname);
  }
  
    
  }
  
}
