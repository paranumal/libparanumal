#include "ins3D.h"

void insRun3D(ins_t *ins, char *options){

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

//   // Set subscycling
//   iint subcycling =0;
//   if(strstr(options,"SUBCYCLING")){ subcycling = 1; }

  occa::initTimer(mesh->device);
  ins->NtimeSteps = 1; // !!!!!!!!!!!!!
  for(iint tstep=0;tstep<ins->NtimeSteps;++tstep){
  #if 0
    // ok it seems 
    if(tstep<100){
      // no advection, first order in time
      ins->b0 = 1.f, ins->a0 = 0.f, ins->c0 = 0.0f; // (1,1,1)
      ins->b1 = 0.f, ins->a1 = 0.f, ins->c1 = 0.0f;
      ins->b2 = 0.f, ins->a2 = 0.f, ins->c2 = 0.0f;
      ins->g0 = 1.f;
    }
    else if(tstep<200){
      // advection, first order in time, no increment
      ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 0.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f;      
    }
    else if(tstep<300){
      // advection, second order in time, first order increment
      ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 0.0f;  // 2
      ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.5f;
    }
    else if(tstep<400)
    {
      // advection, second order in time, first order increment
      ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
      ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.5f;
    }
    else{
      ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 = 1.0f;
      ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 = 0.0f;
      ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
      ins->g0 =  11.f/6.f;
    }
  #else
   if(tstep<1){
       //advection, first order in time, increment
      ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 1.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f;      
    }
    else {
    //advection, second order in time, no increment
    ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
    ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
    ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
    ins->g0 =  1.5f;
     }
    
  #endif

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    
    insAdvectionStep3D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insHelmholtzStep3D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insPoissonStep3D(  ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insUpdateStep3D(   ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);

    printf("tstep = %d\n", tstep);
    if(strstr(options, "REPORT")){
      if(((tstep+1)%(10*ins->errorStep))==0){
        insReport3D(ins, tstep+1,options);
      }
    }
    
   

#if 0 // For time accuracy test fed history with exact solution
    if(tstep<1){
      iint offset = (mesh->Nelements+mesh->totalHaloPairs);
     // Overwrite Velocity
     for(iint e=0;e<mesh->Nelements;++e){
        for(iint n=0;n<mesh->Np;++n){
          const iint id = n + mesh->Np*e;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];

          dfloat a = M_PI/4.0f, d = M_PI/2.0f; 
          dfloat u = -a*( exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y) )* exp(-d*d*t);
          dfloat v = -a*( exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z) )* exp(-d*d*t);
          dfloat w = -a*( exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x) )* exp(-d*d*t);
          dfloat p = -a*a*exp(-2.f*d*d*t)*( exp(2.f*a*x) +exp(2.f*a*y)+exp(2.f*a*z))*( 
                      sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))+
                      sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))+
                      sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))   );   

          // Current time
          const int index0 = (ins->index+0)%3;
          const iint id0   = n + mesh->Np*(e+index0*offset);
          ins->U[id0] = u; 
          ins->V[id0] = v; 
          ins->W[id0] = v; 
          ins->P[id0] = p;
        }
      }
     
       ins->o_U.copyFrom(ins->U);
       ins->o_V.copyFrom(ins->V);
       ins->o_W.copyFrom(ins->W);
       ins->o_P.copyFrom(ins->P);
    }
#endif
  }


 // For Final Time
insReport3D(ins, ins->NtimeSteps+1,options);

#if 0
insErrorNorms2D(ins, ins->finalTime, options);
#endif

//Deallocate Halo MPI storage
free(tSendBuffer);
free(tRecvBuffer);
free(vSendBuffer);
free(vRecvBuffer);
free(pSendBuffer);
free(pRecvBuffer);
}



