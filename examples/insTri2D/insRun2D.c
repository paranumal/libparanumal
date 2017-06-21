#include "ins2D.h"

void insRun2D(ins_t *ins, char *options){

  mesh2D *mesh = ins->mesh;
  // Write Initial Data
  insReport2D(ins, 0, options);
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

  // Set subscycling 
  iint subcycling =0;
  if(strstr(options,"SUBCYCLING")){ subcycling = 1; }
    

  occa::initTimer(mesh->device);

  // for(iint tstep=0;tstep<ins->NtimeSteps;++tstep){
  for(iint tstep=0;tstep<1;++tstep){
  
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
    else if(tstep<400){
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

//     if(tstep<100){
//       // no advection, first order in time
//       ins->b0 = 1.f, ins->a0 = 0.f, ins->c0 = 0.0f; // (1,1,1)
//       ins->b1 = 0.f, ins->a1 = 0.f, ins->c1 = 0.0f;
//       ins->b2 = 0.f, ins->a2 = 0.f, ins->c2 = 0.0f;
//       ins->g0 = 1.f;
//     }
//     #if 0
//     else if(tstep<150){
//       // advection, first order in time, no increment
//       ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 1.0f;  // 2
//       ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
//       ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
//       ins->g0 =  1.f;      
//     }
//     #endif
//     else if(tstep<200){
//       // advection, second order in time, first order increment
//       ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
//       ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
//       ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
//       ins->g0 =  1.5f;
//     }
//     else if(tstep<250){
//       // advection, second order in time, first order increment
//       ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
//       ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
//       ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
//       ins->g0 =  1.5f;
//       }
// #if 1
//     else{
//       // not ok 
//       ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 =  2.0f;
//       ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 =  -1.0f;
//       ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
//       ins->g0 =  11.f/6.f;
//     }
// #endif

  
    // if(tstep<1){
       //advection, first order in time, no increment
      ins->b0 =  1.f,  ins->a0 =  1.0f, ins->c0 = 1.0f;  // 2
      ins->b1 =  0.f,  ins->a1 =  0.0f, ins->c1 = 0.0f; // -1
      ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
      ins->g0 =  1.f;      
     // }
    //   else{ 
    //     // if(tstep<2){
    //   // advection, second order in time, no increment
    //   ins->b0 =  2.f,  ins->a0 =  2.0f, ins->c0 = 1.0f;  // 2
    //   ins->b1 = -0.5f, ins->a1 = -1.0f, ins->c1 = 0.0f; // -1
    //   ins->b2 =  0.f,  ins->a2 =  0.f,  ins->c2 = 0.0f;
    //   ins->g0 =  1.5f;
    // }

    // else{
    //   // 
    //   ins->b0 =  3.f,       ins->a0  =  3.0f, ins->c0 =  1.0f;
    //   ins->b1 = -1.5f,      ins->a1  = -3.0f, ins->c1 =  0.0f;
    //   ins->b2 =  1.f/3.f,   ins->a2  =  1.0f, ins->c2 =  0.0f;
    //   ins->g0 =  11.f/6.f;
    // }


#endif

    ins->lambda = ins->g0 / (ins->dt * ins->nu);

    switch(subcycling){
      case 1:
        insAdvectionSubCycleStep2D(ins, tstep,tSendBuffer,tRecvBuffer,vSendBuffer,vRecvBuffer, options);
      break;

      case 0:
       insAdvectionStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
      break;
    }

    insHelmholtzStep2D(ins, tstep, tHaloBytes,tSendBuffer,tRecvBuffer, options);
    insPoissonStep2D(  ins, tstep, vHaloBytes,vSendBuffer,vRecvBuffer, options);
    insUpdateStep2D(   ins, tstep, pHaloBytes,pSendBuffer,pRecvBuffer, options);

    printf("tstep = %d\n", tstep);
    if(strstr(options, "REPORT")){
      if(((tstep+1)%ins->errorStep)==0){
	      insReport2D(ins, tstep+1,options);
      }
    }
  }


   // For Final Time
  insReport2D(ins, ins->NtimeSteps+1,options);
  
  // Compute the Error and Write To a File
  #if 1
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V);  
  ins->o_P.copyTo(ins->P);

  dfloat t  = ins->dt; 
  dfloat uerr = 0 , verr= 0 , perr = 0; 
  const iint offset = ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);
   for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint id   = n+e*mesh->Np;
      dfloat x  = mesh->x[id];
      dfloat y  = mesh->y[id];
      //
      #if 0
      dfloat uex = -sin(2.0 *M_PI*y)*exp(-ins->nu*4.0*M_PI*M_PI*t);
      dfloat vex =  sin(2.0 *M_PI*x)*exp(-ins->nu*4.0*M_PI*M_PI*t);
      dfloat pex = -cos(2.0 *M_PI*y)*cos(2.0*M_PI*x)*exp(-ins->nu*8.0*M_PI*M_PI*t);
      #else
      dfloat lambda = 1./(2. * ins->nu) - sqrt(1./(4.*ins->nu * ins->nu) + 4.*M_PI*M_PI) ;
      dfloat uex = 1.0 - exp(lambda*x)*cos(2.*M_PI*y);
      dfloat vex =  lambda/(2.*M_PI)*exp(lambda*x)*sin(2.*M_PI*y);
      dfloat pex = 0.5*(1.0- exp(2.*lambda*x));
      #endif


      //
      id += offset;
      uerr = mymax(uerr, fabs(uex - ins->U[id]));
      verr = mymax(verr, fabs(vex - ins->V[id]));
      perr = mymax(perr, fabs(pex - ins->P[id]));

      ins->U[id] = fabs(uex- ins->U[id]);
      ins->V[id] = fabs(vex- ins->V[id]);
      ins->P[id] = fabs(pex- ins->P[id]);
      }
    }

  // compute maximum over all processes
    dfloat guerr, gverr, gperr;
    MPI_Allreduce(&uerr, &guerr, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&verr, &gverr, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&perr, &gperr, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0){
      char fname[BUFSIZ];
      sprintf(fname, "VortexConvergenceTest.txt");
      FILE *fp;
      fp = fopen(fname, "a");
      fprintf(fp,"%d %.5e %.5e %.5e\n", mesh->N, guerr, gverr, gperr);
    }

   
   // Write Error to last file
  ins->o_U.copyFrom(ins->U);
  ins->o_V.copyFrom(ins->V);  
  ins->o_P.copyFrom(ins->P);
  insReport2D(ins, ins->NtimeSteps+2,options);
  #endif


  //





  // Deallocate Halo MPI storage
  free(tSendBuffer);
  free(tRecvBuffer);
  free(vSendBuffer);
  free(vRecvBuffer);
  free(pSendBuffer);
  free(pRecvBuffer);
}



