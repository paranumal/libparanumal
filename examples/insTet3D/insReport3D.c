#include "ins3D.h"

void insReport3D(ins_t *ins, int tstep, char *options){

  dfloat t = (tstep)*ins->dt;
  
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V); 
  ins->o_W.copyTo(ins->W); 
  ins->o_P.copyTo(ins->P);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // do error stuff on host
  insError3D(ins, t, options);

 
  if(strstr(options, "VTU")){ 
    printf("Writing output file\n");
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins3D/");
    // sprintf(fname, "%sfoo_%04d", fname,rank);
    sprintf(fname, "foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    
    insPlotVTU3D(ins, fname);
  } 
}

