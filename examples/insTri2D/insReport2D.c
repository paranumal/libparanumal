#include "ins2D.h"

void insReport2D(ins_t *ins, iint tstep, char *options){

  dfloat t = (tstep+1)*ins->dt;
  
  // copy data back to host
  
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V);  
  ins->o_P.copyTo(ins->P);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // do error stuff on host
  insError2D(ins, t, options);

 
  if(strstr(options, "VTU")){ 
    // compute vorticity
    //insComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
    // output field files
    char fname[BUFSIZ];
    sprintf(fname, "foo_%04d", rank);
    sprintf(fname, "%s_%04d.vtu", fname, tstep/ins->errorStep);
    insPlotVTU2D(ins, fname);
  } 
}

