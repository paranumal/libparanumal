#include "insQuad2D.h"

void insComputeDt2D(ins_t *ins, dfloat currentTime, const char *options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  mesh_t *mesh = ins->mesh;

  dfloat cfl = .2; // pretty good estimate (at least for subcycling LSERK4)

  dfloat dt = 1e9;
  
  // Find Maximum Velocity
  for(dlong e=0;e<mesh->Nelements;++e){

    dfloat hmin = 1e9, hmax = 0;

    for(int f=0;f<mesh->Nfaces;++f){
      dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
      dfloat sJ   = mesh->sgeo[sid + SJID];
      dfloat invJ = mesh->sgeo[sid + IJID];
      
      dfloat hest = 2./(sJ*invJ);
      
      hmin = mymin(hmin, hest);
      hmax = mymax(hmax, hest);
    }

    dfloat umax = .1;
    
    for(int n=0;n<mesh->Np;++n){
      const dlong id = n + mesh->Np*e;

      dfloat uxn = ins->U[id];
      dfloat uyn = ins->V[id];

      //Squared maximum velocity
      dfloat numax = uxn*uxn + uyn*uyn;
      umax = mymax(umax, numax);
    }

    // Maximum Velocity
    umax = sqrt(umax);
    
    dfloat magVel = mymax(umax,1.0); // Correction for initial zero velocity
    
    dfloat localdt = cfl* hmin/( (mesh->N+1.)*(mesh->N+1.) * magVel) ;
    
    dt = mymin(localdt, dt);
  }
  
  // MPI_Allreduce to get global minimum dt
  MPI_Allreduce(&dt, &(ins->dt), 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  if(strstr(options,"SUBCYCLING")){
    ins->dt         = ins->Nsubsteps*ins->dt;
    ins->NtimeSteps = (ins->finalTime-currentTime)/ins->dt;
    ins->dt         = (ins->finalTime-currentTime)/ins->NtimeSteps;
    ins->sdt        = ins->dt/ins->Nsubsteps;
  } else{
    ins->NtimeSteps = (ins->finalTime-currentTime)/ins->dt;
    ins->dt         = (ins->finalTime-currentTime)/ins->NtimeSteps;
  }

  if (rank==0) {
    printf("cfl = %g\n",  cfl);
    printf("dt = %g\n",   ins->dt);
  }
  
  if (strstr(options,"SUBCYCLING")&&rank==0)
    printf("dt: %.8f and sdt: %.8f ratio: %.8f \n", ins->dt, ins->sdt, ins->dt/ins->sdt);

  ins->idt = 1.0/ins->dt;

  if (rank==0)
    printf("Nsteps = %d NerrStep= %d dt = %.8e\n", ins->NtimeSteps,ins->errorStep, ins->dt);

  ins->lambda = ins->g0 / (ins->dt * ins->nu);
  ins->idt = 1.0/ins->dt; 
  ins->ig0 = 1.0/ins->g0; 

  // need to reinterpolate histories ?
  
}
