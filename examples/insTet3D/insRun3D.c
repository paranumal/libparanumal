#include "ins3D.h"

void insRun3D(ins_t *ins, char *options){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mesh3D *mesh = ins->mesh;
  
  // Write Initial Data
  insReport3D(ins, 0, options);
  
  // Allocate MPI buffer for velocity step solver!! May Change Later!!!!!!
  iint tHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
  dfloat *tSendBuffer = (dfloat*) malloc(tHaloBytes);
  dfloat *tRecvBuffer = (dfloat*) malloc(tHaloBytes);

  iint vHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
  dfloat *vSendBuffer = (dfloat*) malloc(vHaloBytes);
  dfloat *vRecvBuffer = (dfloat*) malloc(vHaloBytes);

  // No need to do like this, just for consistency
  iint pHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  dfloat *pSendBuffer = (dfloat*) malloc(pHaloBytes);
  dfloat *pRecvBuffer = (dfloat*) malloc(pHaloBytes);


  occa::initTimer(mesh->device);

  ins->NtimeSteps = 271000; 

  for(iint tstep=0;tstep<ins->NtimeSteps;++tstep){
    if(tstep<1){
      //advection, first order in time, increment
      ins->b0 =  1.f,  ins->a0 =  0.0f, ins->c0 = 1.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f; 
      ins->ExplicitOrder = 1;    
      
      ins->lambda = ins->g0 / (ins->dt * ins->nu);
      ins->idt = 1.0/ins->dt; 
      ins->ig0 = 1.0/ins->g0; 
    // } else if(tstep<2) {
    //   //advection, second order in time, increment
    //   ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
    //   ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
    //   ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
    //   ins->g0 =  1.5f;
    //   ins->ExplicitOrder=2;

    //   ins->lambda = ins->g0 / (ins->dt * ins->nu);
    //   ins->idt = 1.0/ins->dt; 
    //   ins->ig0 = 1.0/ins->g0; 
    // } else {
    //   //advection, third order in time, increment
    //   ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 = 1.0f;
    //   ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 = 0.0f;
    //   ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
    //   ins->g0 =  11.f/6.f;
    //   ins->ExplicitOrder=3;

    //   ins->lambda = ins->g0 / (ins->dt * ins->nu);
    //   ins->idt = 1.0/ins->dt; 
    //   ins->ig0 = 1.0/ins->g0; 
    }

    if(strstr(options,"SUBCYCLING")) {
      insAdvectionSubCycleStep3D(ins, tstep,tSendBuffer,tRecvBuffer,vSendBuffer,vRecvBuffer, options);
    } else {
      insAdvectionStep3D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    }

    insHelmholtzStep3D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insPoissonStep3D(  ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insUpdateStep3D(   ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);

    
    if(((tstep+1)%(ins->errorStep))==0){
      if (rank==0) printf("\rtstep = %d, time = %3.2E, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, (tstep+1)*ins->dt, ins->NiterU, ins->NiterV, ins->NiterW,  ins->NiterP);
      insReport3D(ins, tstep+1,options);
    }

    if (rank==0) printf("\rtstep = %d, time = %3.2E, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, (tstep+1)*ins->dt, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);

     if(strstr(options, "REPORT")){
      if(((tstep+1)%(ins->errorStep))==0){
        if (rank==0)  printf("\rtstep = %d, time = %3.2E, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d \n", tstep+1, (tstep+1)*ins->dt, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
        insErrorNorms3D(ins, (tstep+1)*ins->dt, options);
      }
    }
    
#if 0
    if(tstep<1){
      iint Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
      dfloat tt   = (tstep+1)*ins->dt;
     for(iint e=0;e<mesh->Nelements;++e){
        for(iint n=0;n<mesh->Np;++n){
          iint id = n + mesh->Np*e;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];

          dfloat a = M_PI/4.0f, d = M_PI/2.0f; 
          dfloat u = -a*( exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y) )* exp(-d*d*tt);
          dfloat v = -a*( exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z) )* exp(-d*d*tt);
          dfloat w = -a*( exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x) )* exp(-d*d*tt);
          dfloat p = -a*a*exp(-2.f*d*d*tt)*( exp(2.f*a*x) +exp(2.f*a*y)+exp(2.f*a*z))*( 
                      sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))+
                      sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))+
                      sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))   );   
          // Current time
          id += ins->index*Ntotal; 
          ins->U[id] = u; 
          ins->V[id] = v; 
          ins->W[id] = w; 
          ins->P[id] = p;
        }
      }
     
      ins->o_U.copyFrom(ins->U);
      ins->o_V.copyFrom(ins->V);
      ins->o_W.copyFrom(ins->W);
      ins->o_P.copyFrom(ins->P);
    }
#endif
  }

#if 1
  // For Final Time
  printf("\n");
  insReport3D(ins, ins->NtimeSteps+1,options);
  dfloat finaltime = (ins->NtimeSteps)*ins->dt;
  insErrorNorms3D(ins, ins->finalTime, options);
#endif

  //Deallocate Halo MPI storage
  free(tSendBuffer);
  free(tRecvBuffer);
  free(vSendBuffer);
  free(vRecvBuffer);
  free(pSendBuffer);
  free(pRecvBuffer);
}



