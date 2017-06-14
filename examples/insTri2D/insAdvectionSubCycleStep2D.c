#include "ins2D.h"

// complete a time step using LSERK4
void insAdvectionSubCycleStep2D(ins_t *ins, iint tstep, 
                          dfloat * tSendBuffer, dfloat * tRecvBuffer, 
                          dfloat * sendBuffer, dfloat * recvBuffer, 
                          char   * options){



mesh2D *mesh = ins->mesh; 
// field offset 
iint offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);  

// May be modified : Solver may do smt while exchanging
if(mesh->totalHaloPairs>0){
  ins->totalHaloExtractKernel(mesh->Nelements,
                              mesh->totalHaloPairs,
                              mesh->o_haloElementList,
                              offset,
                              ins->o_U,
                              ins->o_V,
                              ins->o_P,
                              ins->o_tHaloBuffer);
  // copy extracted halo to HOST
  ins->o_tHaloBuffer.copyTo(sendBuffer);
  // start halo exchange
  meshHaloExchangeStart(mesh,
                        mesh->Np*(ins->NTfields)*sizeof(dfloat),
                        sendBuffer,
                        recvBuffer);
}

// Compute Volume Contribution
ins->gradientVolumeKernel(mesh->Nelements,
                          mesh->o_vgeo,
                          mesh->o_DrT,
                          mesh->o_DsT,
                          offset,
                          ins->o_P,
                          ins->o_Px,
                          ins->o_Py);


if(mesh->totalHaloPairs>0){
  meshHaloExchangeFinish(mesh);

  ins->o_tHaloBuffer.copyFrom(recvBuffer);

  ins->totalHaloScatterKernel(mesh->Nelements,
                              mesh->totalHaloPairs,
                              mesh->o_haloElementList,
                              offset,
                              ins->o_U,
                              ins->o_V,
                              ins->o_P,
                              ins->o_tHaloBuffer);

}





// Solve Stokes Problem if Nonlinear solver is deactivated
dfloat activate_advection = 0.f; 
if(ins->a0){activate_advection  = 1.f; }

const iint voffset = 0;   // Velocity halo offset for substepping exchange

// printf("-------------------------------------------------------------------------\n");
iint Ntotal =  (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
ins->o_Ue.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
ins->o_Ve.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
ins->o_Ud.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
ins->o_Vd.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));


for(iint ststep = 0; ststep<ins->Nsubsteps;++ststep){

  dfloat time = tstep*ins->dt + ststep*ins->sdt;    
   // LSERK4 stages
  for(iint rk=0;rk<mesh->Nrk;++rk){
    // intermediate stage time
    dfloat t = time +  ins->sdt*mesh->rkc[rk]; 

    if(mesh->totalHaloPairs>0){
 
    ins->velocityHaloExtractKernel(mesh->Nelements,
                               mesh->totalHaloPairs,
                               mesh->o_haloElementList,
                               voffset, // 0 offset
                               ins->o_Ud,
                               ins->o_Vd,
                               ins->o_vHaloBuffer);

    // copy extracted halo to HOST 
    ins->o_vHaloBuffer.copyTo(sendBuffer);            
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                          mesh->Np*(ins->NVfields)*sizeof(dfloat), 
                          sendBuffer,
                          recvBuffer);
  }

  
    // Compute Volume Contribution
    if(strstr(options, "CUBATURE")){

      ins->subCycleCubatureVolumeKernel(mesh->Nelements,
                 mesh->o_vgeo,
                 mesh->o_cubDrWT,
                 mesh->o_cubDsWT,
                 mesh->o_cubInterpT,
                 ins->o_Ue,
                 ins->o_Ve,
                 ins->o_Ud,
                 ins->o_Vd,
                 ins->o_rhsU,
                 ins->o_rhsV);
    }
    else{
      //Compute Volume
      ins->subCycleVolumeKernel(mesh->Nelements,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                ins->o_Ue,
                                ins->o_Ve,
                                ins->o_Ud,
                                ins->o_Vd,
                                ins->o_rhsU,
                                ins->o_rhsV);

    }


    if(mesh->totalHaloPairs>0){
  
    meshHaloExchangeFinish(mesh);

    ins->o_vHaloBuffer.copyFrom(recvBuffer); 

    ins->velocityHaloScatterKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                voffset, //0 offset
                                ins->o_Ud,
                                ins->o_Vd,
                                ins->o_vHaloBuffer);
  }


  // Compute Volume Contribution
  if(strstr(options, "CUBATURE")){
    ins->subCycleCubatureSurfaceKernel(mesh->Nelements,
                                        mesh->o_sgeo,
                                        mesh->o_intInterpT,
                                        mesh->o_intLIFTT,
                                        mesh->o_vmapM,
                                        mesh->o_vmapP,
                                        mesh->o_EToB,
                                        t,
                                        mesh->o_intx,
                                        mesh->o_inty,
                                        ins->o_Ue,
                                        ins->o_Ve,
                                        ins->o_Ud,
                                        ins->o_Vd,
                                        ins->o_rhsU,
                                        ins->o_rhsV);
    }
  else{
     //Surface Kernel
    ins->subCycleSurfaceKernel(mesh->Nelements,
                              mesh->o_sgeo,
                              mesh->o_LIFTT,
                              mesh->o_vmapM,
                              mesh->o_vmapP,
                              mesh->o_EToB,
                              t,
                              mesh->o_x,
                              mesh->o_y,
                              ins->o_Ue,
                              ins->o_Ve,
                              ins->o_Ud,
                              ins->o_Vd,
                              ins->o_rhsU,
                              ins->o_rhsV);

  }

  // Update Kernel
  ins->subCycleRKUpdateKernel(mesh->Nelements,
                              activate_advection,
                              ins->sdt,
                              mesh->rka[rk],
                              mesh->rkb[rk],
                              ins->o_rhsU,
                              ins->o_rhsV,
                              ins->o_resU, 
                              ins->o_resV,
                              ins->o_Ud,
                              ins->o_Vd);
  }


//printf("Extrapolating Velocity to %d \n", ststep+1);
// Extrapolate Velocity
iint offset1 = mesh->Nelements+mesh->totalHaloPairs;
ins->subCycleExtKernel((mesh->Nelements+mesh->totalHaloPairs),
                      ststep,
                      ins->sdt,
                      ins->dt,
                      ins->index,
                      offset1,
                      ins->o_U,
                      ins->o_V,
                      ins->o_Ue,
                      ins->o_Ve);


 }



  //copy into next time level storage 
  iint index1 = (ins->index+1)%3;
  ins->o_Ud.copyTo(ins->o_NU,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
  ins->o_Vd.copyTo(ins->o_NV,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);


 dfloat t = tstep*ins->dt;  
// Compute Surface Conribution
  ins->gradientSurfaceKernel(mesh->Nelements,
           mesh->o_sgeo,
           mesh->o_LIFTT,
           mesh->o_vmapM,
           mesh->o_vmapP,
           mesh->o_EToB,
           mesh->o_x,
           mesh->o_y,
           t,
           ins->dt,
           ins->a0,
           ins->a1,
           ins->a2,
           ins->index,
           mesh->Nelements+mesh->totalHaloPairs,
           0, 
           ins->o_PI, //not used
           ins->o_P,
           ins->o_Px,
           ins->o_Py);  




}