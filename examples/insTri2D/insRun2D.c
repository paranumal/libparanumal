#include "ins2D.h"



void insRun2D(ins_t *ins, char *options){

  mesh2D *mesh = ins->mesh;
  
  // Write Initial Data
  insReport2D(ins, 0, options);
  // Allocate MPI buffer for velocity step solver
  iint  tHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
  dfloat  *tSendBuffer = (dfloat*) malloc(tHaloBytes);
  dfloat  *tRecvBuffer = (dfloat*) malloc(tHaloBytes);

  iint vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
  dfloat *vSendBuffer = (dfloat*) malloc(vHaloBytes);
  dfloat *vRecvBuffer = (dfloat*) malloc(vHaloBytes);

  // No need to do like this, just for consistency
  iint pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  dfloat *pSendBuffer = (dfloat*) malloc(pHaloBytes);
  dfloat *pRecvBuffer = (dfloat*) malloc(pHaloBytes);

  // Set subscycling
  iint subcycling =0;
  if(strstr(options,"SUBCYCLING")){ subcycling = 1; }

  occa::initTimer(mesh->device);

  occaTimerTic(mesh->device,"INS");
  
  // if(ins->Nsubsteps)
  // ins->NtimeSteps = 160/ins->Nsubsteps;
  // else
  ins->NtimeSteps=100;
  
  double tic_tot = 0.f, elp_tot = 0.f; 
  double tic_adv = 0.f, elp_adv = 0.f;
  double tic_pre = 0.f, elp_pre = 0.f;
  double tic_vel = 0.f, elp_vel = 0.f;
  double tic_upd = 0.f, elp_upd = 0.f;

  // MPI_Barrier(MPI_COMM_WORLD); 
  tic_tot = MPI_Wtime(); 
  for(iint tstep=0;tstep<ins->NtimeSteps;++tstep){
    if(tstep<1){
       //advection, first order in time, increment
      ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 1.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f; 
      ins->ExplicitOrder = 1;    
      
      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    } else { //if(tstep<2) 
      //advection, second order in time, no increment
      ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
      ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.5f;
      ins->ExplicitOrder=2;

      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    }
    // else{
    // //advection, second order in time, no increment
    // ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 = 1.0f;
    // ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 = 0.0f;
    // ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
    // ins->g0 =  11.f/6.f;
    // ins->ExplicitOrder=3;
    // }

   
    tic_adv = MPI_Wtime(); 
    // if(strstr(options,"ALGEBRAIC")){
     switch(subcycling){
      case 1:
        occaTimerTic(mesh->device,"AdvectionSubStep");
        insAdvectionSubCycleStep2D(ins, tstep,tSendBuffer,tRecvBuffer,vSendBuffer,vRecvBuffer, options);
        occaTimerToc(mesh->device,"AdvectionSubStep");
      break;

      case 0:
      occaTimerTic(mesh->device,"AdvectionStep");
      insAdvectionStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
      occaTimerToc(mesh->device,"AdvectionStep");
      break;
    }

    elp_adv += (MPI_Wtime() - tic_adv); 


    tic_vel = MPI_Wtime();  

    occaTimerTic(mesh->device,"HelmholtzStep");
    insHelmholtzStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options); 
    occaTimerToc(mesh->device,"HelmholtzStep");

    elp_vel += (MPI_Wtime() - tic_vel);

    tic_pre = MPI_Wtime();  

    occaTimerTic(mesh->device,"PoissonStep");
    insPoissonStep2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    occaTimerToc(mesh->device,"PoissonStep");

    elp_pre += (MPI_Wtime() - tic_pre);


    tic_upd = MPI_Wtime();
    occaTimerTic(mesh->device,"UpdateStep");
    insUpdateStep2D(ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);
    occaTimerToc(mesh->device,"UpdateStep"); 

    elp_upd += (MPI_Wtime() - tic_upd);
    
    // } 

    // else {
    // insAdvectionStepSS2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    // insPoissonStepSS2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    // insHelmholtzStepSS2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    // }

    // printf("tstep = %d of %d\n", tstep,ins->NtimeSteps);
    
    #if 1
    occaTimerTic(mesh->device,"Report");
    if(strstr(options, "VTU")){
      if(((tstep+1)%(ins->errorStep))==0){
        printf("\n");
        insReport2D(ins, tstep+1,options);
      }
    }

    printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);

     if(strstr(options, "REPORT")){
      if(((tstep+1)%(ins->errorStep))==0){
        printf("\n");
        insErrorNorms2D(ins, (tstep+1)*ins->dt, options);
      }
    }
    
    occaTimerToc(mesh->device,"Report");
    #endif

    
    #if 0// For time accuracy test fed history with exact solution
        if(tstep<1){
          iint Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
          dfloat tt = (tstep+1)*ins->dt;
         // Overwrite Velocity
         for(iint e=0;e<mesh->Nelements;++e){
            for(iint n=0;n<mesh->Np;++n){
              iint id = n + mesh->Np*e;
              dfloat x = mesh->x[id];
              dfloat y = mesh->y[id];

              dfloat u0 = -sin(2.0 *M_PI*y)*exp(-ins->nu*4.0*M_PI*M_PI*tt); ;
              dfloat v0 =  sin(2.0 *M_PI*x)*exp(-ins->nu*4.0*M_PI*M_PI*tt); 
              dfloat p0 = -cos(2.0 *M_PI*y)*cos(2.f*M_PI*x)*exp(-ins->nu*8.f*M_PI*M_PI*tt);

              id += ins->index*Ntotal; 

              ins->U[id] = u0; 
              ins->V[id] = v0; 
              ins->P[id] = p0;
            }
          }
           ins->o_U.copyFrom(ins->U);
           ins->o_V.copyFrom(ins->V);
           ins->o_P.copyFrom(ins->P);
        }
    #endif
  }

occaTimerToc(mesh->device,"INS");

// MPI_Barrier(MPI_COMM_WORLD);
elp_tot = MPI_Wtime() - tic_tot; 




 // compute maximum over all processes
  double gelp_tot  = 0.f, gelp_adv = 0.f, gelp_vel = 0.f , gelp_pre = 0.f , gelp_upd = 0.f ;
 
  MPI_Allreduce(&elp_tot, &gelp_tot, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elp_adv, &gelp_adv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elp_vel, &gelp_vel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elp_pre, &gelp_pre, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&elp_upd, &gelp_upd, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);



int rank, size;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

if(rank==0){
  printf("%2d %2d %.5e %.5e %.5e %.5e %.5e\n", mesh->N,size,gelp_tot, gelp_adv, gelp_vel, gelp_pre, gelp_upd); 
  
  char fname[BUFSIZ]; sprintf(fname, "insScaling2D.dat");
  FILE *fp; fp = fopen(fname, "a");
  fprintf(fp, "%2d %2d %.5e %.5e %.5e %.5e %.5e\n", mesh->N,size,gelp_tot, gelp_adv, gelp_vel, gelp_pre, gelp_upd); 
  fclose(fp);
}




// if(rank==0)
// printf("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", gelp_tot/double(size), gelp_adv/double(size), 
//                                          gelp_vel/double(size), gelp_pre/double(size), gelp_upd/double(size));



#if 1
dfloat finalTime = ins->NtimeSteps*ins->dt;
insReport2D(ins, ins->NtimeSteps,options);
insErrorNorms2D(ins, finalTime, options);
#endif

// Deallocate Halo MPI storage
free(tSendBuffer);
free(tRecvBuffer);
free(vSendBuffer);
free(vRecvBuffer);
free(pSendBuffer);
free(pRecvBuffer);




  //int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0) occa::printTimer();

}



