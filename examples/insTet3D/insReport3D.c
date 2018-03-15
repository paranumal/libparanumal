#include "ins3D.h"

void insReport3D(ins_t *ins, int tstep, char *options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh3D *mesh = ins->mesh;

  dfloat t = (tstep)*ins->dt;
  
  int offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  ins->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_DrT,
                       mesh->o_DsT,
                       mesh->o_DtT,
                       offset,
                       ins->o_U,
                       ins->o_V,
                       ins->o_W,
                       ins->o_Vx,
                       ins->o_Vy,
                       ins->o_Vz);

  ins->divergenceKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_DrT,
                       mesh->o_DsT,
                       mesh->o_DtT,
                       offset,
                       ins->o_U,
                       ins->o_V,
                       ins->o_W,
                       ins->o_Div);

  // gather-scatter
  ellipticParallelGatherScatterTet3D(mesh, mesh->ogs, ins->o_Vx, dfloatString, "add");
  ellipticParallelGatherScatterTet3D(mesh, mesh->ogs, ins->o_Vy, dfloatString, "add");  
  ellipticParallelGatherScatterTet3D(mesh, mesh->ogs, ins->o_Vz, dfloatString, "add");  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vx, ins->o_Vx);
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vy, ins->o_Vy);
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vz, ins->o_Vz);  

  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V); 
  ins->o_W.copyTo(ins->W); 
  ins->o_P.copyTo(ins->P);

  ins->o_Vx.copyTo(ins->Vx);
  ins->o_Vy.copyTo(ins->Vy);
  ins->o_Vz.copyTo(ins->Vz);
  ins->o_Div.copyTo(ins->Div);
  
  insError3D(ins, t, options);
  
  if (rank==0) printf("Writing output file\n");
  
  // output field files
  char fname[BUFSIZ];

  if(strstr(options, "VTU")){   
    // sprintf(fname, "/u0/outputs/ins3D/");
    // sprintf(fname, "%sfoo_%04d", fname,rank);
    sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    
    insPlotVTU3D(ins, fname);
  } else if(strstr(options, "SLICE")){   
    //slice data (cylinders)
    // const int Nslices = 4;
    // const char *sliceDim[4] = {"x","y","y","z"};
    // const dfloat sliceX[4] = {0.0,1.0,-1.0,20.};

    //slice data (channel)
    // const int Nslices = 4;
    // const char *sliceDim[4] = {"x","x","y","z"};
    // const dfloat sliceX[4] = {0.001,5,0.0,0.0};

    //slice data (fence)
    // const int Nslices = 4;
    // const char *sliceDim[4] = {"x","x","y","z"};
    // const dfloat sliceX[4] = {0.0,5,0.25,0.0};

    //slice data (cub)
    const int Nslices = 3;
    const char *sliceDim[4] = {"x","y","z"};
    const dfloat sliceX[4] = {0.0,0.0,0.0};

    // output field files
    sprintf(fname, "slice_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    //sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    insPlotSlice3D(ins, fname,Nslices, sliceDim,sliceX);
  } else if(strstr(options, "CONTOUR")){ 
  
    sprintf(fname, "contour_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    //sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank,tstep/ins->errorStep);
    insPlotContour3D(ins, fname, options);
  } 
}

