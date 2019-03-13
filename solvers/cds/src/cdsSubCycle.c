/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "cds.h"

// complete a time step using LSERK4
void cdsSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_U, occa::memory o_S, occa::memory o_Sd){


   mesh_t *mesh = cds->mesh;
  const dlong NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);  

  //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0){
#if 0
    cds->haloExtractKernel(mesh->Nelements,
                           mesh->totalHaloPairs,
                           mesh->o_haloElementList,
                           cds->vOffset,
                           cds->sOffset,
                           o_U,
                           o_S,
                           cds->o_haloBuffer);

    // copy extracted halo to HOST 
    cds->o_haloBuffer.copyTo(cds->sendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
        mesh->Np*(cds->NSfields+cds->NVfields)*sizeof(dfloat),
        cds->sendBuffer,
        cds->recvBuffer);

    meshHaloExchangeFinish(mesh);

    cds->o_haloBuffer.copyFrom(cds->recvBuffer); 

    cds->haloScatterKernel(mesh->Nelements,
         mesh->totalHaloPairs,
         cds->vOffset,
         cds->sOffset,
         o_U,
         o_S,
         cds->o_haloBuffer);
  #else
     cds->haloGetKernel(mesh->totalHaloPairs,
                         cds->vOffset,
                         cds->sOffset,
                         mesh->o_haloElementList,
                         mesh->o_haloGetNodeIds,
                         o_U,
                         o_S,
                         cds->o_haloBuffer);
          
      dlong Ndata = (cds->NVfields+cds->NSfields)*mesh->Nfp*mesh->totalHaloPairs;
          // copy extracted halo to HOST 
      cds->o_haloBuffer.copyTo(cds->sendBuffer, Ndata*sizeof(dfloat), 0);// zero offset

       // start halo exchange
      meshHaloExchangeStart(mesh,
                            mesh->Nfp*(cds->NVfields+ cds->NSfields)*sizeof(dfloat),
                            cds->sendBuffer,
                            cds->recvBuffer);

      meshHaloExchangeFinish(mesh);

      cds->o_haloBuffer.copyFrom(cds->recvBuffer, Ndata*sizeof(dfloat), 0);  

      cds->haloPutKernel(mesh->totalHaloPairs,
                          cds->vOffset,
                          cds->sOffset,
                          mesh->o_haloElementList,
                          mesh->o_haloPutNodeIds,
                          o_U,
                          o_S,
                          cds->o_haloBuffer);
  #endif
  }


  const dfloat tn0 = time - 0*cds->dt;
  const dfloat tn1 = time - 1*cds->dt;
  const dfloat tn2 = time - 2*cds->dt;


  dfloat zero = 0.0, one = 1.0;
  int izero = 0;
  dfloat b, bScale=0;

   // Solve for Each SubProblem
  for (int torder=(cds->ExplicitOrder-1); torder>=0; torder--){
    
    b=cds->extbdfB[torder];
    bScale += b; 

    // printf("Torder = %d Writing bScale: %.4f \n", torder, bScale);

    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder*cds->NSfields*cds->Ntotal;

    if (torder==cds->ExplicitOrder-1) { //first substep
      cds->scaledAddKernel(cds->NSfields*cds->Ntotal, b, toffset, o_S, zero, izero, o_Sd);
    } else { //add the next field
      cds->scaledAddKernel(cds->NSfields*cds->Ntotal, b, toffset, o_S,  one, izero, o_Sd);
    }     

    // SubProblem  starts from here from t^(n-torder)
    const dfloat tsub = time - torder*cds->dt;
    // Advance SubProblem to t^(n-torder+1) 
    for(int ststep = 0; ststep<cds->Nsubsteps;++ststep){
     
      const dfloat tstage = tsub + ststep*cds->sdt;     
     
      for(int rk=0;rk<cds->SNrk;++rk){// LSERK4 stages
        // Extrapolate velocity to subProblem stage time
        dfloat t = tstage +  cds->sdt*cds->Srkc[rk]; 

        switch(cds->ExplicitOrder){
          case 1:
            cds->extC[0] = 1.f; cds->extC[1] = 0.f; cds->extC[2] = 0.f;
            break;
          case 2:
            cds->extC[0] = (t-tn1)/(tn0-tn1);
            cds->extC[1] = (t-tn0)/(tn1-tn0);
            cds->extC[2] = 0.f; 
            break;
          case 3:
            cds->extC[0] = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
            cds->extC[1] = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
            cds->extC[2] = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
            break;
        }
        cds->o_extC.copyFrom(cds->extC);

        //compute advective velocity fields at time t
        cds->subCycleExtKernel(NtotalElements,
                               Nstages,
                               cds->vOffset,
                               cds->o_extC,
                               o_U,
                               cds->o_Ue);

       if(mesh->totalHaloPairs>0){
          // make sure compute device is ready to perform halo extract
          mesh->device.finish();
          // switch to data stream
          mesh->device.setStream(mesh->dataStream);

#if 0          

          cds->scalarHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 cds->sOffset, 
                                 o_Sd,
                                 cds->o_shaloBuffer);

          // copy extracted halo to HOST 
          cds->o_shaloBuffer.copyTo(cds->ssendBuffer,"async: true");   
#else
      cds->scalarHaloGetKernel(mesh->totalHaloPairs,
                                   cds->sOffset,
                                   mesh->o_haloElementList,
                                   mesh->o_haloGetNodeIds,
                                   o_S,
                                   cds->o_shaloBuffer);
          
      dlong Ndata = (cds->NSfields)*mesh->Nfp*mesh->totalHaloPairs;
          // copy extracted halo to HOST 
      cds->o_shaloBuffer.copyTo(cds->ssendBuffer, Ndata*sizeof(dfloat), 0, "async: true");
#endif


          mesh->device.setStream(mesh->defaultStream);
        }

        // Compute Volume Contribution
        occaTimerTic(mesh->device,"AdvectionVolume");        
        if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          cds->subCycleCubatureVolumeKernel(mesh->Nelements,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDWmatrices,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              cds->vOffset,
              cds->sOffset,             
              cds->o_Ue,
                   o_Sd,
              cds->o_rhsSd);
        } else{
          cds->subCycleVolumeKernel(mesh->Nelements,
                                    mesh->o_vgeo,
                                    mesh->o_Dmatrices,
                                    cds->vOffset,
                                    cds->sOffset,           
                                    cds->o_Ue,
                                    o_Sd,
                                    cds->o_rhsSd);

        }
        occaTimerToc(mesh->device,"AdvectionVolume");

        if(mesh->totalHaloPairs>0){
          // make sure compute device is ready to perform halo extract
          mesh->device.setStream(mesh->dataStream);
          mesh->device.finish();

#if 0

          // start halo exchange
          meshHaloExchangeStart(mesh,
                                mesh->Np*(cds->NSfields)*sizeof(dfloat), 
                                cds->ssendBuffer,
                                cds->srecvBuffer);
        

          meshHaloExchangeFinish(mesh);

          cds->o_shaloBuffer.copyFrom(cds->srecvBuffer,"async: true"); 

          cds->scalarHaloScatterKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 cds->sOffset,
                                 o_Sd,
                                 cds->o_shaloBuffer);
#else
            // start halo exchange
      meshHaloExchangeStart(mesh,
                            mesh->Nfp*(cds->NSfields)*sizeof(dfloat),
                            cds->ssendBuffer,
                            cds->srecvBuffer);

      meshHaloExchangeFinish(mesh);

      dlong Ndata = (cds->NSfields)*mesh->Nfp*mesh->totalHaloPairs;
      cds->o_shaloBuffer.copyFrom(cds->srecvBuffer, Ndata*sizeof(dfloat), 0, "async: true");  

      cds->scalarHaloPutKernel(mesh->totalHaloPairs,
                          cds->sOffset,
                          mesh->o_haloElementList,
                          mesh->o_haloPutNodeIds,
                          o_S,
                          cds->o_shaloBuffer);

#endif
          mesh->device.finish();
          
          mesh->device.setStream(mesh->defaultStream);
          mesh->device.finish();
        }

        //Surface Kernel
        occaTimerTic(mesh->device,"AdvectionSurface");
        if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          cds->subCycleCubatureSurfaceKernel(mesh->Nelements,
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
               cds->vOffset,
               cds->sOffset,               
               cds->o_Ue,
                    o_Sd,
               cds->o_rhsSd);
        } else{
          cds->subCycleSurfaceKernel(mesh->Nelements,
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
             cds->vOffset,
             cds->sOffset,             
             cds->o_Ue,
                  o_Sd,
             cds->o_rhsSd);
        }
        occaTimerToc(mesh->device,"AdvectionSurface");
          
        // Update Kernel
        occaTimerTic(mesh->device,"AdvectionUpdate");
        cds->subCycleRKUpdateKernel(mesh->Nelements,
                                    cds->sdt,
                                    cds->Srka[rk],
                                    cds->Srkb[rk],
                                    cds->sOffset,
                                    cds->o_rhsSd,
                                    cds->o_resS, 
                                         o_Sd);
        occaTimerToc(mesh->device,"AdvectionUpdate");
      }
    }
  }

}

