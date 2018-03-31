#include "insHex3D.h"

// complete a time step using LSERK4
void insAdvectionSubCycleStepHex3D(ins_t *ins, int tstep, char   * options){
 
  //printf("SUBSTEP METHOD : SEMI-LAGRAGIAN OIFS METHOD\n");
  mesh3D *mesh = ins->mesh;

  const dlong NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);
  const dlong Ntotal         =  NtotalElements*mesh->Np;  

  const dlong voffset = 0; 
  // field offset at this step
  dlong offset0 = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  
  //Exctract Halo On Device, Assumes History is already done!
  if(mesh->totalHaloPairs>0){
    ins->totalHaloExtractKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                offset0,
                                ins->o_U,
                                ins->o_V,
                                ins->o_W,
                                ins->o_P,
                                ins->o_tHaloBuffer);

    // copy extracted halo to HOST
    ins->o_tHaloBuffer.copyTo(ins->tSendBuffer);

    // start halo exchange    
    meshHaloExchangeStart(mesh,
                          mesh->Np*(ins->NTfields)*sizeof(dfloat),
                          ins->tSendBuffer,
                          ins->tRecvBuffer);
  }

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    ins->o_tHaloBuffer.copyFrom(ins->tRecvBuffer);
    ins->totalHaloScatterKernel(mesh->Nelements,
                                mesh->totalHaloPairs,
                                mesh->o_haloElementList,
                                offset0,
                                ins->o_U,
                                ins->o_V,
                                ins->o_W,
                                ins->o_P,
                                ins->o_tHaloBuffer);
  }



  if (ins->a0) {// skip if nonlinear term is deactivated (Stokes solves)
    const dfloat tn0 = (tstep-0)*ins->dt;
    const dfloat tn1 = (tstep-1)*ins->dt;
    const dfloat tn2 = (tstep-2)*ins->dt;

    //storage for subcycled velocity fields
    // use NU and NV 
    occa::memory o_Ud = ins->o_NU;
    occa::memory o_Vd = ins->o_NV;
    occa::memory o_Wd = ins->o_NW;

    // construct interpolating lagrange polynomial
    dfloat c0 = 0.f, c1 = 0.f, c2 = 0.f;
    
    dfloat zero = 0.0, one = 1.0;
    int izero = 0;

    dfloat b, bScale=0;

    // Solve for Each SubProblem
    for (int torder=ins->ExplicitOrder-1; torder>=0; torder--){
      
      if (torder==2) b=ins->b2;
      if (torder==1) b=ins->b1;
      if (torder==0) b=ins->b0;

      bScale += b;

      // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
      const int sindex = (ins->index + 3 - torder)%3; 
      dlong offset = sindex*Ntotal;
      if (torder==ins->ExplicitOrder-1) { //first substep
        ins->scaledAddKernel(mesh->Nelements*mesh->Np, b, offset, ins->o_U, zero, izero, o_Ud);
        ins->scaledAddKernel(mesh->Nelements*mesh->Np, b, offset, ins->o_V, zero, izero, o_Vd);
        ins->scaledAddKernel(mesh->Nelements*mesh->Np, b, offset, ins->o_W, zero, izero, o_Wd);
      } else { //add the next field
        ins->scaledAddKernel(mesh->Nelements*mesh->Np, b, offset, ins->o_U,  one, izero, o_Ud);
        ins->scaledAddKernel(mesh->Nelements*mesh->Np, b, offset, ins->o_V,  one, izero, o_Vd);
        ins->scaledAddKernel(mesh->Nelements*mesh->Np, b, offset, ins->o_W,  one, izero, o_Wd);
      }     

      // SubProblem  starts from here from t^(n-torder)
      const dfloat tsub = tstep*ins->dt - torder*ins->dt;
      // Advance SubProblem to t^(n-torder+1) 
      for(int ststep = 0; ststep<ins->Nsubsteps;++ststep){
        const dfloat tstage = tsub + ststep*ins->sdt;     
        for(int rk=0;rk<mesh->Nrk;++rk){// LSERK4 stages
          // Extrapolate velocity to subProblem stage time
          dfloat t = tstage +  ins->sdt*mesh->rkc[rk]; 

          switch(ins->ExplicitOrder){
            case 1:
              c0 = 1.f; c1 = 0.f; c2 = 0.f;
              break;
            case 2:
              c0 = (t-tn1)/(tn0-tn1);
              c1 = (t-tn0)/(tn1-tn0);
              c2 = 0.f; 
              break;
            case 3:
              c0 = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
              c1 = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
              c2 = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
              break;
          }

          //compute advective velocity fields at time t
          ins->subCycleExtKernel(NtotalElements,
                                 ins->index,
                                 NtotalElements,
                                 c0,
                                 c1,
                                 c2,
                                 ins->o_U,
                                 ins->o_V,
                                 ins->o_W,
                                 ins->o_Ue,
                                 ins->o_Ve,
                                 ins->o_We);

          if(mesh->totalHaloPairs>0){
            // make sure compute device is ready to perform halo extract
            mesh->device.finish();

            // switch to data stream
            mesh->device.setStream(mesh->dataStream);

            ins->velocityHaloExtractKernel(mesh->Nelements,
                                     mesh->totalHaloPairs,
                                     mesh->o_haloElementList,
                                     voffset, // 0 offset
                                     o_Ud,
                                     o_Vd,
                                     o_Wd,
                                     ins->o_vHaloBuffer);

            // copy extracted halo to HOST 
            ins->o_vHaloBuffer.asyncCopyTo(ins->vSendBuffer);            
            mesh->device.setStream(mesh->defaultStream);        
          }
          occaTimerTic(mesh->device,"AdvectionVolume");
          
          // Compute Volume Contribution
          if(strstr(options, "CUBATURE")){
            ins->subCycleCubatureVolumeKernel(mesh->Nelements,
                                              mesh->o_vgeo,
                                              mesh->o_cubvgeo,
                                              mesh->o_cubDWT,
                                              mesh->o_cubInterpT,
                                              mesh->o_cubProjectT,
                                              ins->o_Ue,
                                              ins->o_Ve,
                                              ins->o_We,
                                                   o_Ud,
                                                   o_Vd,
                                                   o_Wd,
                                              ins->o_cU,
                                              ins->o_cV,
                                              ins->o_cW,
                                              ins->o_cUd,
                                              ins->o_cVd,
                                              ins->o_cWd,
                                              ins->o_rhsU,
                                              ins->o_rhsV,
                                              ins->o_rhsW);
          } else{
            ins->subCycleVolumeKernel(mesh->Nelements,
                                      mesh->o_vgeo,
                                      mesh->o_D,
                                      ins->o_Ue,
                                      ins->o_Ve,
                                      ins->o_We,
                                           o_Ud,
                                           o_Vd,
                                           o_Wd,
                                      ins->o_rhsU,
                                      ins->o_rhsV,
                                      ins->o_rhsW);

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

            ins->o_vHaloBuffer.asyncCopyFrom(ins->vRecvBuffer);  

            ins->velocityHaloScatterKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      mesh->o_haloElementList,
                                      voffset, //0 offset
                                         o_Ud,
                                         o_Vd,
                                         o_Wd,
                                      ins->o_vHaloBuffer);
            mesh->device.finish();
            
            mesh->device.setStream(mesh->defaultStream);
            mesh->device.finish();
          }

          //Surface Kernel
          occaTimerTic(mesh->device,"AdvectionSurface");
          if(strstr(options, "CUBATURE")){
            ins->subCycleCubatureSurfaceKernel(mesh->Nelements,
                                                mesh->o_vgeo,
                                                mesh->o_cubsgeo,
                                                mesh->o_vmapM,
                                                mesh->o_vmapP,
                                                mesh->o_EToB,
                                                mesh->o_cubInterpT,
                                                mesh->o_cubProjectT,
                                                bScale,
                                                t,
                                                mesh->o_intx,
                                                mesh->o_inty,
                                                mesh->o_intz,
                                                ins->o_Ue,
                                                ins->o_Ve,
                                                ins->o_We,
                                                     o_Ud,
                                                     o_Vd,
                                                     o_Wd,
                                                ins->o_rhsU,
                                                ins->o_rhsV,
                                                ins->o_rhsW);
          } else{
            ins->subCycleSurfaceKernel(mesh->Nelements,
                                      mesh->o_sgeo,
                                      mesh->o_vmapM,
                                      mesh->o_vmapP,
                                      mesh->o_EToB,
                                      bScale,
                                      t,
                                      mesh->o_x,
                                      mesh->o_y,
                                      mesh->o_z,
                                      ins->o_Ue,
                                      ins->o_Ve,
                                      ins->o_We,
                                           o_Ud,
                                           o_Vd,
                                           o_Wd,
                                      ins->o_rhsU,
                                      ins->o_rhsV,
                                      ins->o_rhsW);
          }
          occaTimerToc(mesh->device,"AdvectionSurface");
            
          // Update Kernel
          occaTimerTic(mesh->device,"AdvectionUpdate");
          ins->subCycleRKUpdateKernel(mesh->Nelements,
                                ins->sdt,
                                mesh->rka[rk],
                                mesh->rkb[rk],
                                ins->o_rhsU,
                                ins->o_rhsV,
                                ins->o_rhsW,
                                ins->o_resU, 
                                ins->o_resV,
                                ins->o_resW,
                                     o_Ud,
                                     o_Vd,
                                     o_Wd);
          occaTimerToc(mesh->device,"AdvectionUpdate");
        }
      }
    }
  }

  dfloat tp1 = (tstep+1)*ins->dt;

  // Compute Volume Contribution for Pressure
  occaTimerTic(mesh->device,"GradientVolume");
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_D,
                            offset0,
                            ins->o_P,
                            ins->o_Px,
                            ins->o_Py,
                            ins->o_Pz);
  occaTimerToc(mesh->device,"GradientVolume");
 
  if (strstr(ins->pSolverOptions,"IPDG")) {
    const int solverid = 0; // Pressure Solve
    occaTimerTic(mesh->device,"GradientSurface");
    // Compute Surface Conribution
    ins->gradientSurfaceKernel(mesh->Nelements,
                               mesh->o_sgeo,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               mesh->o_EToB,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               tp1,
                               ins->dt,
                               ins->c0,
                               ins->c1,
                               ins->c2,
                               ins->index,
                               mesh->Nelements+mesh->totalHaloPairs,
                               solverid, // pressure BCs
                               ins->o_PI, //not used
                               ins->o_P,
                               ins->o_Px,
                               ins->o_Py,
                               ins->o_Pz);
    occaTimerToc(mesh->device,"GradientSurface");
  }
}