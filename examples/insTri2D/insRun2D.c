#include "ins2D.h"

void insRun2D(solver_t *ins, char *options){

  mesh2D *mesh = ins->mesh; 

  // Allocate MPI send buffer
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*ins->Nfields*sizeof(dfloat);
  dfloat *sendBuffer = (dfloat*) malloc(haloBytes);
  dfloat *recvBuffer = (dfloat*) malloc(haloBytes);

  occa::initTimer(mesh->device);

 
  //First Order Coefficients
  ins->a0 = 1.0, ins->b0 = 1.0,  ins->a1 = 0.0, ins->b1 =0.0, ins->g0 = 1.0; 

   
  insAdvectionStep2D(ins, 0, haloBytes, sendBuffer, recvBuffer,options);




  
  // // Switch to second order
  // ins->a0 = 2.0, ins->b0 = 2.0,  ins->a1 = -0.5, ins->b1 =-1.0, ins->g0 = 1.5; 


  // occa::initTimer(mesh->device);

  // for(iint tstep=1;tstep<ins->NtimeSteps;++tstep){
  //      //
  //      insAdvectionStep2D(ins, tstep, haloBytes, sendBuffer, recvBuffer,options);
  //      //
  //      //insPressureStep2D(ins, tstep, haloBytes, sendBuffer, recvBuffer,options);









  //    if(strstr(options, "REPORT")){
  //     if((tstep%ins->errorStep)==0){
  //       insReport2D(ins, tstep,options);
  //     }
  //   }
  //   }
  
  // // For Final Time
  insReport2D(ins, ins->NtimeSteps,options);

  occa::printTimer();

  // Deallocate Halo MPI storage
  free(recvBuffer);
  free(sendBuffer);
}



