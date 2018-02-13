#include "ins3D.h"

void insReport3D(ins_t *ins, iint tstep, char *options){

  dfloat t = (tstep)*ins->dt;
  
  if(strstr(options, "VTU")){ 
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
 
    if (rank==0) printf("Writing output file\n");
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins3D/");
    // sprintf(fname, "%sfoo_%04d", fname,rank);
    sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    
    insPlotVTU3D(ins, fname);
  } else if(strstr(options, "SLICE")){ 
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
 
    if (rank==0) printf("Writing output file\n");
    
    //slice data
    const int Nslices = 4;
    const char *sliceDim[4] = {"x","y","y","z"};
    const dfloat sliceX[4] = {0.0,1.0,-1.0,20.};

    // output field files
    char fname[BUFSIZ];
    sprintf(fname, "slice_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    //sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    insPlotSlice3D(ins, fname,Nslices, sliceDim,sliceX);
  } if(strstr(options, "ADAPTIVECONTOUR")){ 
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
 
    if (rank==0) printf("Writing output file\n");
    
    // output field files
    char fname[BUFSIZ];
    sprintf(fname, "/scratch/contour_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    //sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    insPlotAdaptiveContour3D(ins, fname);
  } 
}

