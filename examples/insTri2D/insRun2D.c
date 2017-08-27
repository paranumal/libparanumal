#include "ins2D.h"

void insRun2D(ins_t *ins, char *options){

  mesh2D *mesh = ins->mesh;
  
  // Write Initial Data
  // insReport2D(ins, 0, options);
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



 
 // occa::initTimer(mesh->device);
 ins->NtimeSteps = 10;
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
      ins->ExplicitOrder = 1;     
   }
    else
     // if(tstep<2) 
    {
    //advection, second order in time, no increment
    ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
    ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
    ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
    ins->g0 =  1.5f;
    ins->ExplicitOrder=2;
    }
    // else{
    // //advection, second order in time, no increment
    // ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 = 1.0f;
    // ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 = 0.0f;
    // ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
    // ins->g0 =  11.f/6.f;
    // ins->ExplicitOrder=3;
    // }


    
  #endif
    ins->lambda = ins->g0 / (ins->dt * ins->nu);

    if(strstr(options,"ALGEBRAIC")){
     switch(subcycling){
      case 1:
        insAdvectionSubCycleStep2D(ins, tstep,tSendBuffer,tRecvBuffer,vSendBuffer,vRecvBuffer, options);
      break;

      case 0:
      // VOLUME KERNELS
      // mesh->device.finish();
      // occa::tic("AdvectionStep");
      printf("Without substepping.....\n");
      insAdvectionStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
      // mesh->device.finish();
      // occa::toc("AdvectionStep");  
      break;
    }

    // mesh->device.finish();
    // occa::tic("HelmholtzStep");
    insHelmholtzStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    // mesh->device.finish();
    // occa::toc("HelmholtzStep");  

    // mesh->device.finish();
    // occa::tic("PoissonStep");
    insPoissonStep2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    // mesh->device.finish();
    // occa::toc("PoissonStep");  
    
    // mesh->device.finish();
    // occa::tic("UpdateStep");
    insUpdateStep2D(ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);
    // mesh->device.finish();
    // occa::toc("UpdateStep");  
    } else {

    insAdvectionStepSS2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insPoissonStepSS2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insHelmholtzStepSS2D(ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);





    }

    printf("tstep = %d of %d\n", tstep,ins->NtimeSteps);

    if(strstr(options, "REPORT")){
      if(((tstep+1)%(ins->errorStep))==0){
        insReport2D(ins, tstep+1,options);
        insErrorNorms2D(ins, (tstep+1)*ins->dt, options);
      }
    }
    
#if 1 // For time accuracy test fed history with exact solution
    if(tstep<1){
      iint Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
      dfloat tt = (tstep+1)*ins->dt;
      dfloat *tmpU = (dfloat*) calloc(Ntotal, sizeof(dfloat));
      dfloat *tmpV = (dfloat*) calloc(Ntotal, sizeof(dfloat));
      dfloat *tmpP = (dfloat*) calloc(Ntotal, sizeof(dfloat));
     // Overwrite Velocity
     for(iint e=0;e<mesh->Nelements;++e){
        for(iint n=0;n<mesh->Np;++n){
          const iint id = n + mesh->Np*e;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];

          dfloat u0 = -sin(2.0 *M_PI*y)*exp(-ins->nu*4.0*M_PI*M_PI*tt); ;
          dfloat v0 =  sin(2.0 *M_PI*x)*exp(-ins->nu*4.0*M_PI*M_PI*tt); 
          dfloat p0 = -cos(2.0 *M_PI*y)*cos(2.f*M_PI*x)*exp(-ins->nu*8.f*M_PI*M_PI*tt);

          tmpU[id] = u0; 
          tmpV[id] = v0; 
          tmpP[id] = p0;
        }
      }
       ins->o_U.copyFrom(tmpU,Ntotal*sizeof(dfloat),ins->index*Ntotal*sizeof(dfloat));
       ins->o_V.copyFrom(tmpV,Ntotal*sizeof(dfloat),ins->index*Ntotal*sizeof(dfloat));
       ins->o_P.copyFrom(tmpP,Ntotal*sizeof(dfloat),Ntotal*sizeof(dfloat));

       free(tmpU);
       free(tmpV);
       free(tmpP);       
    }
#endif



  }


 // For Final Time
insReport2D(ins, ins->NtimeSteps,options);

#if 0
dfloat finalTime = ins->NtimeSteps*ins->dt;
insErrorNorms2D(ins, finalTime, options);
#endif

// Deallocate Halo MPI storage
free(tSendBuffer);
free(tRecvBuffer);
free(vSendBuffer);
free(vRecvBuffer);
free(pSendBuffer);
free(pRecvBuffer);

// occa::printTimer();

}



