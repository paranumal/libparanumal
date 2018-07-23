#include "ins.h"

// complete a time step using LSERK4
void insSubCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_U, occa::memory o_Ud){
 
  //printf("SUBSTEP METHOD : SEMI-LAGRAGIAN OIFS METHOD\n");
  mesh_t *mesh = ins->mesh;

  const dlong NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);  

  //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0){
    ins->velocityHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ins->fieldOffset,
                                 o_U,
                                 ins->o_vHaloBuffer);

    // copy extracted halo to HOST 
    ins->o_vHaloBuffer.copyTo(ins->vSendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*(ins->NVfields)*sizeof(dfloat),
                         ins->vSendBuffer,
                         ins->vRecvBuffer);

    meshHaloExchangeFinish(mesh);

    ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

    ins->velocityHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  ins->fieldOffset,
                                  o_U,
                                  ins->o_vHaloBuffer);
  }

  
  const dfloat tn0 = time - 0*ins->dt;
  const dfloat tn1 = time - 1*ins->dt;
  const dfloat tn2 = time - 2*ins->dt;

  dfloat zero = 0.0, one = 1.0;
  int izero = 0;

  dfloat b, bScale=0;

  // Solve for Each SubProblem
  for (int torder=ins->ExplicitOrder-1; torder>=0; torder--){
    
    b=ins->extbdfB[torder];
    bScale += b;

    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder*ins->NVfields*ins->Ntotal;

    if (torder==ins->ExplicitOrder-1) { //first substep
      ins->scaledAddKernel(ins->NVfields*ins->Ntotal, b, toffset, o_U, zero, izero, o_Ud);
    } else { //add the next field
      ins->scaledAddKernel(ins->NVfields*ins->Ntotal, b, toffset, o_U,  one, izero, o_Ud);
    }     

    // SubProblem  starts from here from t^(n-torder)
    const dfloat tsub = time - torder*ins->dt;
    // Advance SubProblem to t^(n-torder+1) 
    for(int ststep = 0; ststep<ins->Nsubsteps;++ststep){
      const dfloat tstage = tsub + ststep*ins->sdt;     
      for(int rk=0;rk<mesh->Nrk;++rk){// LSERK4 stages
        // Extrapolate velocity to subProblem stage time
        dfloat t = tstage +  ins->sdt*mesh->rkc[rk]; 

        switch(ins->ExplicitOrder){
          case 1:
            ins->extC[0] = 1.f; ins->extC[1] = 0.f; ins->extC[2] = 0.f;
            break;
          case 2:
            ins->extC[0] = (t-tn1)/(tn0-tn1);
            ins->extC[1] = (t-tn0)/(tn1-tn0);
            ins->extC[2] = 0.f; 
            break;
          case 3:
            ins->extC[0] = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
            ins->extC[1] = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
            ins->extC[2] = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
            break;
        }
        ins->o_extC.copyFrom(ins->extC);

        //compute advective velocity fields at time t
        ins->subCycleExtKernel(NtotalElements,
                               Nstages,
                               ins->fieldOffset,
                               ins->o_extC,
                               o_U,
                               ins->o_Ue);

        if(mesh->totalHaloPairs>0){
          // make sure compute device is ready to perform halo extract
          mesh->device.finish();

          // switch to data stream
          mesh->device.setStream(mesh->dataStream);

          ins->velocityHaloExtractKernel(mesh->Nelements,
                                   mesh->totalHaloPairs,
                                   mesh->o_haloElementList,
                                   ins->fieldOffset, 
                                   o_Ud,
                                   ins->o_vHaloBuffer);

          // copy extracted halo to HOST 
          ins->o_vHaloBuffer.copyTo(ins->vSendBuffer,"async: true");            
          mesh->device.setStream(mesh->defaultStream);
        }

        // Compute Volume Contribution
        occaTimerTic(mesh->device,"AdvectionVolume");        
        if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          ins->subCycleCubatureVolumeKernel(mesh->Nelements,
                     mesh->o_vgeo,
                     mesh->o_cubvgeo,
                     mesh->o_cubDWmatrices,
                     mesh->o_cubInterpT,
                     mesh->o_cubProjectT,
                     ins->fieldOffset,
                     ins->o_Ue,
                          o_Ud,
                     ins->o_cU,     
                     ins->o_cUd,     
                     ins->o_rhsUd);
        } else{
          ins->subCycleVolumeKernel(mesh->Nelements,
                                    mesh->o_vgeo,
                                    mesh->o_Dmatrices,
                                    ins->fieldOffset,
                                    ins->o_Ue,
                                         o_Ud,
                                    ins->o_rhsUd);

        }
        occaTimerToc(mesh->device,"AdvectionVolume");

        if(mesh->totalHaloPairs>0){
          // make sure compute device is ready to perform halo extract
          mesh->device.setStream(mesh->dataStream);
          mesh->device.finish();

          // start halo exchange
          meshHaloExchangeStart(mesh,
                              mesh->Np*(ins->NVfields)*sizeof(dfloat), 
                              ins->vSendBuffer,
                              ins->vRecvBuffer);
        

          meshHaloExchangeFinish(mesh);

          ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer,"async: true"); 

          ins->velocityHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    ins->fieldOffset, //0 ins->fieldOffset
                                    o_Ud,
                                    ins->o_vHaloBuffer);
          mesh->device.finish();
          
          mesh->device.setStream(mesh->defaultStream);
          mesh->device.finish();
        }

        //Surface Kernel
        occaTimerTic(mesh->device,"AdvectionSurface");
        if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          ins->subCycleCubatureSurfaceKernel(mesh->Nelements,
                                              mesh->o_vgeo,
                                              mesh->o_sgeo,
                                              mesh->o_cubsgeo,
                                              mesh->o_intInterpT,
                                              mesh->o_intLIFTT,
                                              mesh->o_cubInterpT,
                                              mesh->o_cubProjectT,
                                              mesh->o_vmapM,
                                              mesh->o_vmapP,
                                              mesh->o_EToB,
                                              bScale,
                                              t,
                                              mesh->o_intx,
                                              mesh->o_inty,
                                              mesh->o_intz,
                                              ins->fieldOffset,
                                              ins->o_Ue,
                                                   o_Ud,
                                              ins->o_rhsUd);
        } else{
          ins->subCycleSurfaceKernel(mesh->Nelements,
                                    mesh->o_sgeo,
                                    mesh->o_LIFTT,
                                    mesh->o_vmapM,
                                    mesh->o_vmapP,
                                    mesh->o_EToB,
                                    bScale,
                                    t,
                                    mesh->o_x,
                                    mesh->o_y,
                                    mesh->o_z,
                                    ins->fieldOffset,
                                    ins->o_Ue,
                                         o_Ud,
                                    ins->o_rhsUd);
        }
        occaTimerToc(mesh->device,"AdvectionSurface");
          
        // Update Kernel
        occaTimerTic(mesh->device,"AdvectionUpdate");
        ins->subCycleRKUpdateKernel(mesh->Nelements,
                              ins->sdt,
                              mesh->rka[rk],
                              mesh->rkb[rk],
                              ins->fieldOffset,
                              ins->o_rhsUd,
                              ins->o_resU, 
                                   o_Ud);
        occaTimerToc(mesh->device,"AdvectionUpdate");
      }
    }
  }
}