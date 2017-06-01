#include "ins2D.h"

void insRun2D(ins_t *ins, char *options){

  mesh2D *mesh = ins->mesh; 
 
 #if 0
  
  solver_t *solver = ins->prsolver;

  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  
  // load rhs into r
  dfloat *cf = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(iint n=0;n<mesh->Np;++n){
      dfloat xn = mesh->x[n+e*mesh->Np];
      dfloat yn = mesh->y[n+e*mesh->Np];
      nrhs[n] = -(2*M_PI*M_PI+1.0)*cos(M_PI*xn)*cos(M_PI*yn);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->Np;++m){
        rhs += mesh->MM[n+m*mesh->Np]*nrhs[m];
      }
      iint id = n+e*mesh->Np;
      
      r[id] = -rhs*J;
      x[id] = 0;
      ins->Pr[id] = nrhs[n];
    }
  }
  free(nrhs);
  free(cf);

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  ellipticSolveTri2D(solver,1.0 , o_r, o_x, ins->prsolverOptions);

  // copy solution from DEVICE to HOST
  o_x.copyTo(ins->Pr);

  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat exact = cos(M_PI*xn)*cos(M_PI*yn);
      dfloat error = fabs(exact-ins->Pr[id]);
      
      maxError = mymax(maxError, error);
    }
  }
  
  dfloat globalMaxError = maxError;
  //MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  //if(rank==0)
  printf("globalMaxError = %g\n", globalMaxError);
  
  //meshPlotVTU2D(mesh, "foo", 0);

#endif
 insReport2D(ins, 0, options);
  // Allocate MPI buffer for velocity step solver
  iint helmholtzHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NTfields)*sizeof(dfloat);
  dfloat *helmholtzSendBuffer = (dfloat*) malloc(helmholtzHaloBytes);
  dfloat *helmholtzRecvBuffer = (dfloat*) malloc(helmholtzHaloBytes);
  
  //
  iint poissonHaloBytes = mesh->totalHaloPairs*mesh->Np*(ins->NVfields)*sizeof(dfloat);
  dfloat *poissonSendBuffer = (dfloat*) malloc(poissonHaloBytes);
  dfloat *poissonRecvBuffer = (dfloat*) malloc(poissonHaloBytes);

  // No need to do like this, just for consistency 
  iint updateHaloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  dfloat *updateSendBuffer = (dfloat*) malloc(updateHaloBytes);
  dfloat *updateRecvBuffer = (dfloat*) malloc(updateHaloBytes);

  occa::initTimer(mesh->device);
  //First Order Coefficients
  ins->a0 = 1.0, ins->b0 = 1.f;
  ins->a1 = 0.f, ins->b1 = 0.f;
  ins->b2 = 0.f, ins->b2 = 0.f;
  ins->g0 = 1.f; 
  
  ins->lamda = ins->g0 / (ins->dt * ins->nu); // Update Lamda for first order
  // One First Order Step 
  insHelmholtzStep2D(ins, 0, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options);
  //insPoissonStep2D(ins,   0, poissonHaloBytes  , poissonSendBuffer  , poissonRecvBuffer  , options);
  //insUpdateStep2D(ins, 0, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);
  //
 //  printf("Finished first order step\n");
  
 //  // Switch to second order
 //  ins->a0 =  2.0,  ins->b0  = 2.0f;
 //  ins->a1 = -0.5f, ins->b1 = -1.0f;
 //  ins->b2 =  0.f,  ins->b2  = 0.f;
 //  ins->g0 =  1.5f; 

 //  ins->lamda = ins->g0 / (ins->dt * ins->nu);
 // // One First Order Step 
 //  insHelmholtzStep2D(ins, 1, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options);
 //  insPoissonStep2D(ins,   1, poissonHaloBytes  , poissonSendBuffer  , poissonRecvBuffer  , options);
 //  insUpdateStep2D(ins, 1, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);
 //  //
 //  printf("Finished second order step\n");



 //  // Switch to third order
 //  ins->a0 =  3.f,       ins->b0  =  3.0f;
 //  ins->a1 = -1.5f,      ins->b1  = -3.0f;
 //  ins->b2 =  1.f/3.f,   ins->b2  =  1.0f;
 //  ins->g0 =  11.f/6.f; 

 //  // // Set-up First Order Solver
 //  ins->lamda = ins->g0 / (ins->dt * ins->nu);
  
 //    for(iint tstep=2;tstep<ins->NtimeSteps;++tstep){
 //      //
 //      insHelmholtzStep2D(ins, tstep, helmholtzHaloBytes, helmholtzSendBuffer, helmholtzRecvBuffer, options);
 //      insPoissonStep2D(ins,   tstep, poissonHaloBytes  , poissonSendBuffer  , poissonRecvBuffer  , options);
 //      insUpdateStep2D(ins, tstep, updateHaloBytes, updateSendBuffer, updateRecvBuffer, options);
 //     if(strstr(options, "REPORT")){
 //      if((tstep%ins->errorStep)==0){
 //        insReport2D(ins, tstep,options);
 //       }
 //     }
 //  }
  
 //  // // For Final Time
 //  insReport2D(ins, ins->NtimeSteps,options);

 //  occa::printTimer();

  // // Deallocate Halo MPI storage
  free(helmholtzSendBuffer);
  free(helmholtzRecvBuffer);
  // //
  free(poissonSendBuffer);
  free(poissonRecvBuffer);
  // //
  free(updateSendBuffer);
  free(updateRecvBuffer);
  //

}



