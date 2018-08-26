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

#include "bns.h"
#define BNS_ANSYC 0

void bnsMRSAABStep(bns_t *bns, int tstep, int haloBytes,
       dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options){


mesh_t *mesh = bns->mesh; 

const dlong offset    = mesh->Np*mesh->Nelements*bns->Nfields;
const dlong pmloffset = mesh->Np*mesh->pmlNelements*bns->Nfields;

  for (int Ntick=0; Ntick < pow(2,mesh->MRABNlevels-1);Ntick++) {

    // intermediate stage time
    dfloat t = bns->dt*(tstep*pow(2.0,mesh->MRABNlevels-1) + Ntick);

    int mrab_order = 0; 

    if(tstep==0)  mrab_order = 0; // first order
    else if(tstep==1) mrab_order = 1; // second order
    else mrab_order = 2; // third order 

    // COMPUTE RAMP FUNCTION
    dfloat fx, fy, fz, intfx, intfy, intfz;
    bnsBodyForce(t, &fx, &fy, &fz, &intfx, &intfy, &intfz);

    int lev;
    for (lev=0;lev<mesh->MRABNlevels;lev++)
      if (Ntick % (1<<lev) != 0) break; //find the max lev to compute rhs
    
      if(mesh->totalHaloPairs>0){
#if BNS_ASYNC 
        mesh->device.setStream(dataStream);


        int Nentries = mesh->Nfp*bns->Nfields*mesh->Nfaces;
        mesh->haloExtractKernel(mesh->totalHaloPairs,
                                Nentries,
                                mesh->o_haloElementList,
                                bns->o_fQM,
                                mesh->o_haloBuffer);

        // copy extracted halo to HOST
        mesh->o_haloBuffer.copyTo(sendBuffer,"async: true");
        mesh->device.setStream(defaultStream);
#else
        int Nentries = mesh->Nfp*bns->Nfields*mesh->Nfaces;
        mesh->haloExtractKernel(mesh->totalHaloPairs,
                                Nentries,
                                mesh->o_haloElementList,
                                bns->o_fQM,
                                mesh->o_haloBuffer);

        mesh->o_haloBuffer.copyTo(sendBuffer);
        // start halo exchange
        meshHaloExchangeStart(mesh, mesh->Nfields*mesh->Nfp*mesh->Nfaces*sizeof(dfloat),
                              sendBuffer, recvBuffer);
#endif
      }


      occaTimerTic(mesh->device, "VolumeKernel");  

      for (int l=0;l<lev;l++) {

        if (mesh->MRABNelements[l]){
          occaTimerTic(mesh->device, "NonPmlVolumeKernel"); 
          bns->volumeKernel(mesh->MRABNelements[l],
                            mesh->o_MRABelementIds[l],
                            offset,
                            mesh->MRABshiftIndex[l],
                            fx, fy, fz,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            bns->o_q,
                            bns->o_rhsq);
                       
          occaTimerToc(mesh->device, "NonPmlVolumeKernel"); 
        }

        if (mesh->MRABpmlNelements[l]){
          occaTimerTic(mesh->device, "PmlVolumeKernel");

          if(bns->pmlcubature){
            bns->pmlVolumeKernel(mesh->MRABpmlNelements[l],
                                 mesh->o_MRABpmlElementIds[l],
                                 mesh->o_MRABpmlIds[l],
                                 offset,   
                                 pmloffset, 
                                 mesh->MRABshiftIndex[l],
                                 fx, fy, fz,
                                 mesh->o_vgeo,
                                 mesh->o_Dmatrices,
                                 bns->o_q,
                                 bns->o_pmlqx,
                                 bns->o_pmlqy,
                                 bns->o_pmlqz,
                                 bns->o_rhsq,
                                 bns->o_pmlrhsqx,
                                 bns->o_pmlrhsqy,
                                 bns->o_pmlrhsqz);
        }else{
            bns->pmlVolumeKernel(mesh->MRABpmlNelements[l],
                                 mesh->o_MRABpmlElementIds[l],
                                 mesh->o_MRABpmlIds[l],
                                 offset,   
                                 pmloffset, 
                                 mesh->MRABshiftIndex[l],
                                 fx, fy, fz,
                                 mesh->o_vgeo,
                                 bns->o_pmlSigmaX,
                                 bns->o_pmlSigmaY,
                                 bns->o_pmlSigmaZ,
                                 mesh->o_Dmatrices,
                                 bns->o_q,
                                 bns->o_pmlqx,
                                 bns->o_pmlqy,
                                 bns->o_pmlqz,
                                 bns->o_rhsq,
                                 bns->o_pmlrhsqx,
                                 bns->o_pmlrhsqy,
                                 bns->o_pmlrhsqz);
        }
      occaTimerToc(mesh->device, "PmlVolumeKernel");
      }
    }

    occaTimerToc(mesh->device, "VolumeKernel");   

     
    occaTimerTic(mesh->device, "RelaxationKernel");
    for (int l=0;l<lev;l++) {
      if (mesh->MRABNelements[l]){
        occaTimerTic(mesh->device,"NonPmlRelaxationKernel");
        bns->relaxationKernel(mesh->MRABNelements[l],
                              mesh->o_MRABelementIds[l],
                              mesh->o_vgeo,
                              mesh->o_cubvgeo,
                              offset,
                              mesh->MRABshiftIndex[l],
                              mesh->o_cubInterpT,
                              mesh->o_cubProjectT,
                              bns->o_q,
                              bns->o_rhsq); 
        occaTimerToc(mesh->device,"NonPmlRelaxationKernel");
      } 

      if (mesh->MRABpmlNelements[l]){
        occaTimerTic(mesh->device,"PmlRelaxationKernel");
        if(bns->pmlcubature){
          bns->pmlRelaxationKernel(mesh->MRABpmlNelements[l],
                                  mesh->o_MRABpmlElementIds[l],
                                  mesh->o_MRABpmlIds[l],
                                  mesh->o_vgeo,
                                  mesh->o_cubvgeo,
                                  offset,   
                                  pmloffset, 
                                  mesh->MRABshiftIndex[l],
                                  mesh->o_cubInterpT,
                                  mesh->o_cubProjectT,
                                  bns->o_pmlSigmaX,
                                  bns->o_pmlSigmaY,
                                  bns->o_pmlSigmaZ,
                                  bns->o_q,
                                  bns->o_pmlqx,
                                  bns->o_pmlqy,
                                  bns->o_pmlqz,
                                  bns->o_rhsq,
                                  bns->o_pmlrhsqx,
                                  bns->o_pmlrhsqy,
                                  bns->o_pmlrhsqz);
        }else{
          bns->pmlRelaxationKernel(mesh->MRABpmlNelements[l],
                                  mesh->o_MRABpmlElementIds[l],
                                  mesh->o_MRABpmlIds[l],
                                  mesh->o_vgeo,
                                  mesh->o_cubvgeo,
                                  offset,   
                                  pmloffset, 
                                  mesh->MRABshiftIndex[l],
                                  mesh->o_cubInterpT,
                                  mesh->o_cubProjectT,
                                  bns->o_q,
                                  bns->o_rhsq);
        }
        occaTimerToc(mesh->device,"PmlRelaxationKernel");
      }
    }
    occaTimerToc(mesh->device, "RelaxationKernel");


    if(mesh->totalHaloPairs>0){
#if BNS_ASYNC 
        mesh->device.setStream(dataStream);
        //make sure the async copy is finished
        mesh->device.finish();

        // start halo exchange
        meshHaloExchangeStart(mesh,
                              mesh->Nfields*mesh->Nfp*mesh->Nfaces*sizeof(dfloat),
                              sendBuffer,
                              recvBuffer);

        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);

        // copy halo data to DEVICE
        size_t foffset = mesh->Nfaces*mesh->Nfp*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
        bns->o_fQM.copyFrom(recvBuffer, haloBytes, foffset,"async: true");
        mesh->device.finish();
        mesh->device.setStream(defaultStream);
#else
        // wait for halo data to arrive
        meshHaloExchangeFinish(mesh);

        // copy halo data to DEVICE
        size_t foffset = mesh->Nfaces*mesh->Nfp*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
        bns->o_fQM.copyFrom(recvBuffer, haloBytes, foffset);

#endif
    }


    // SURFACE KERNELS for boltzmann Nodal DG
    for (int l=0;l<lev;l++) {
      occaTimerTic(mesh->device,"SurfaceKernel");
      if (mesh->MRABNelements[l]){
        occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
        bns->surfaceKernel(mesh->MRABNelements[l],
                            mesh->o_MRABelementIds[l],
                            offset,
                            mesh->MRABshiftIndex[l],
                            t,
                            intfx, intfy, intfz,
                            mesh->o_sgeo,
                            mesh->o_LIFTT,
                            mesh->o_vmapM,
                            mesh->o_mapP,
                            mesh->o_EToB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            bns->o_q,
                            bns->o_fQM,
                            bns->o_rhsq);
        occaTimerToc(mesh->device,"NonPmlSurfaceKernel");
      }

      if (mesh->MRABpmlNelements[l]){
        occaTimerTic(mesh->device,"PmlSurfaceKernel");
        bns->pmlSurfaceKernel(mesh->MRABpmlNelements[l],
                              mesh->o_MRABpmlElementIds[l],
                              mesh->o_MRABpmlIds[l],
                              offset,
                              pmloffset,
                              mesh->MRABshiftIndex[l],
                              t,
                              intfx, intfy, intfz,
                              mesh->o_sgeo,
                              mesh->o_LIFTT,
                              mesh->o_vmapM,
                              mesh->o_mapP,
                              mesh->o_EToB,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              bns->o_q,
                              bns->o_fQM,
                              bns->o_rhsq,
                              bns->o_pmlrhsqx,
                              bns->o_pmlrhsqy,
                              bns->o_pmlrhsqz);
        occaTimerToc(mesh->device,"PmlSurfaceKernel");
      }
      occaTimerToc(mesh->device,"SurfaceKernel");
    }


    for (lev=0;lev<mesh->MRABNlevels;lev++)
      if ((Ntick+1) % (1<<lev) !=0) break; //find the max lev to update
      for (int l=0; l<lev; l++) {

        const int id = mrab_order*mesh->MRABNlevels*3 + l*3;
        occaTimerTic(mesh->device,"UpdateKernel");

        if (mesh->MRABNelements[l]){
          occaTimerTic(mesh->device,"NonPmlUpdateKernel");
            bns->updateKernel(mesh->MRABNelements[l],
                              mesh->o_MRABelementIds[l],
                              offset,
                              mesh->MRABshiftIndex[l],
                              bns->MRSAAB_C[l], //
                              bns->MRAB_A[id+0], //
                              bns->MRAB_A[id+1],
                              bns->MRAB_A[id+2], //
                              bns->MRSAAB_A[id+0], //
                              bns->MRSAAB_A[id+1],
                              bns->MRSAAB_A[id+2], //
                              mesh->o_vmapM,
                              bns->o_rhsq,
                              bns->o_fQM,
                              bns->o_q);
          occaTimerToc(mesh->device,"NonPmlUpdateKernel");
        }

        if (mesh->MRABpmlNelements[l]){
          occaTimerTic(mesh->device,"PmlUpdateKernel");
          bns->pmlUpdateKernel(mesh->MRABpmlNelements[l],
                              mesh->o_MRABpmlElementIds[l],
                              mesh->o_MRABpmlIds[l],
                              offset,
                              pmloffset,
                              mesh->MRABshiftIndex[l],
                              bns->MRSAAB_C[l], //
                              bns->MRAB_A[id+0], //
                              bns->MRAB_A[id+1],
                              bns->MRAB_A[id+2], //
                              bns->MRSAAB_A[id+0], //
                              bns->MRSAAB_A[id+1],
                              bns->MRSAAB_A[id+2], //
                              mesh->o_vmapM,
                              bns->o_rhsq,
                              bns->o_pmlrhsqx,
                              bns->o_pmlrhsqy,
                              bns->o_pmlrhsqz,
                              bns->o_q,
                              bns->o_pmlqx,
                              bns->o_pmlqy,
                              bns->o_pmlqz,
                              bns->o_fQM);
          occaTimerToc(mesh->device,"PmlUpdateKernel");
        }

        occaTimerToc(mesh->device,"UpdateKernel");
        //rotate index
        mesh->MRABshiftIndex[l] = (mesh->MRABshiftIndex[l]+1)%3;
      }


      if (lev<mesh->MRABNlevels) {    
      // const int id = mrab_order*mesh->MRABNlevels*3 + (lev-1)*3; // !!!!!
        const int id = mrab_order*mesh->MRABNlevels*3 + (lev)*3;
        occaTimerTic(mesh->device,"TraceUpdateKernel");

        if (mesh->MRABNhaloElements[lev]){
          occaTimerTic(mesh->device,"NonPmlTraceUpdateKernel");
          bns->traceUpdateKernel(mesh->MRABNhaloElements[lev],
                                mesh->o_MRABhaloIds[lev],
                                offset,
                                mesh->MRABshiftIndex[lev],
                                bns->MRSAAB_C[lev-1], //
                                bns->MRAB_B[id+0], //
                                bns->MRAB_B[id+1],
                                bns->MRAB_B[id+2], //
                                bns->MRSAAB_B[id+0], //
                                bns->MRSAAB_B[id+1],
                                bns->MRSAAB_B[id+2], 
                                mesh->o_vmapM,
                                bns->o_q,
                                bns->o_rhsq,
                                bns->o_fQM);
          occaTimerToc(mesh->device,"NonPmlTraceUpdateKernel");
      }

      if (mesh->MRABpmlNhaloElements[lev]){
        occaTimerTic(mesh->device,"PmlTraceUpdateKernel");
        bns->traceUpdateKernel(mesh->MRABpmlNhaloElements[lev],
                               mesh->o_MRABpmlHaloElementIds[lev],
                               offset,
                               mesh->MRABshiftIndex[lev],
                               bns->MRSAAB_C[lev-1], //
                               bns->MRAB_B[id+0], //
                               bns->MRAB_B[id+1],
                               bns->MRAB_B[id+2], //
                               bns->MRSAAB_B[id+0], //
                               bns->MRSAAB_B[id+1],
                               bns->MRSAAB_B[id+2], 
                               mesh->o_vmapM,
                               bns->o_q,
                               bns->o_rhsq,
                               bns->o_fQM);
        occaTimerToc(mesh->device,"PmlTraceUpdateKernel");
      }
      occaTimerToc(mesh->device,"TraceUpdateKernel");
    }
  }

  
}
