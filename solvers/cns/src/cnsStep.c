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

#include "cns.h"

#define USE_OLD_HALO 0

void cnsDopriStep(cns_t *cns, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = cns->mesh;
  
  //RK step
  for(int rk=0;rk<cns->Nrk;++rk){


    mesh->device.setStream(mesh->defaultStream);
    
    // t_rk = t + C_rk*dt
    dfloat currentTime = time + cns->rkC[rk]*mesh->dt;

    dfloat fx, fy, fz, intfx, intfy, intfz;
    cnsBodyForce(currentTime , &fx, &fy, &fz, &intfx, &intfy, &intfz);

    //compute RK stage 
    // rkq = q + dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
    cns->rkStageKernel(mesh->Nelements,
                       rk,
                       mesh->dt,
                       cns->o_rkA,
                       cns->o_q,
                       cns->o_rkrhsq,
                       cns->o_rkq);
    
    //compute RHS
    // rhsq = F(currentTIme, rkq)

    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nfields;

#if (USE_OLD_HALO)

      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_rkq, cns->o_haloBuffer);
      
      // copy extracted halo to HOST 
      cns->o_haloBuffer.copyTo(cns->sendBuffer);      
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
#else

      // make sure rkq is updated
      mesh->device.finish();  
      
      // launch haloExtractKernel on 2nd stream
      mesh->device.setStream(mesh->dataStream);         

      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_rkq, cns->o_haloBuffer);

      // launch async copy on 2nd stream
      cns->o_haloBuffer.copyTo(cns->sendBuffer,"async: true");

      // switch to stream0 for Volume kernel
      mesh->device.setStream(mesh->defaultStream);      
#endif
    }

    // now compute viscous stresses
    cns->stressesVolumeKernel(mesh->Nelements, 
                              mesh->o_vgeo, 
                              mesh->o_Dmatrices,
                              cns->mu,
                              cns->o_rkq, 
                              cns->o_viscousStresses);

    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){

#if (!USE_OLD_HALO)
      // switch dev->currentStream
      mesh->device.setStream(mesh->dataStream);

      // flush 2nd stream, ensure send buffer is loaded
      mesh->device.finish();

      // run remaining work on stream0
      mesh->device.setStream(mesh->defaultStream);  

      // start halo exchange on default stream
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
#endif
      
      meshHaloExchangeFinish(mesh);
      
      // offset for halo daat
      size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); 
      cns->o_rkq.copyFrom(cns->recvBuffer, cns->haloBytes, offset);

    }

    cns->stressesSurfaceKernel(mesh->Nelements, 
                               mesh->o_sgeo, 
                               mesh->o_LIFTT,
                               mesh->o_vmapM, 
                               mesh->o_vmapP, 
                               mesh->o_EToB, 
                               currentTime,
                               mesh->o_x, 
                               mesh->o_y,
                               mesh->o_z, 
                               cns->mu,
			       intfx, intfy, intfz,
                               cns->o_rkq, 
                               cns->o_viscousStresses);

    // extract stresses halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nstresses;

#if (USE_OLD_HALO)
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
      
      // copy extracted halo to HOST 
      cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
      
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
#else
      // launch haloExtractKernel on 2nd stream
      mesh->device.setStream(mesh->dataStream);                       

      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);

      // launch async copy on 2nd stream
      cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer,"async: true");
      
      // switch to stream0 for Volume kernel
      mesh->device.setStream(mesh->defaultStream);                   
#endif
				 

      
    }

    // compute volume contribution to DG cns RHS
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {
      cns->cubatureVolumeKernel(mesh->Nelements, 
                                cns->advSwitch,
				fx, fy, fz,
                                mesh->o_vgeo,
                                mesh->o_cubvgeo, 
                                mesh->o_cubDWmatrices,
                                mesh->o_cubInterpT,
                                mesh->o_cubProjectT,
                                cns->o_viscousStresses, 
                                cns->o_rkq, 
                                cns->o_rhsq);
    } else {
      cns->volumeKernel(mesh->Nelements, 
                        cns->advSwitch,
			fx, fy, fz,
                        mesh->o_vgeo, 
                        mesh->o_Dmatrices,
                        cns->o_viscousStresses, 
                        cns->o_rkq, 
                        cns->o_rhsq);
    }

    // wait for halo stresses data to arrive
    if(mesh->totalHaloPairs>0){
#if (!USE_OLD_HALO)
      // switch dev->currentStream
      mesh->device.setStream(mesh->dataStream);

      // flush 2nd stream, ensure send buffer is loaded
      mesh->device.finish();

      // run remaining work on stream0
      mesh->device.setStream(mesh->defaultStream);  
      
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
#endif

      meshHaloExchangeFinish(mesh);
	    
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
      
    }

    // compute surface contribution to DG cns RHS (LIFTT ?)
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {
      cns->cubatureSurfaceKernel(mesh->Nelements, 
                                 cns->advSwitch,
                                 mesh->o_vgeo, 
                                 mesh->o_cubsgeo, 
                                 mesh->o_vmapM, 
                                 mesh->o_vmapP, 
                                 mesh->o_EToB,
                                 mesh->o_intInterpT,
                                 mesh->o_intLIFTT, 
                                 currentTime, 
                                 mesh->o_intx, 
                                 mesh->o_inty,
                                 mesh->o_intz, 
                                 cns->mu,
				 intfx, intfy, intfz,
                                 cns->o_rkq, 
                                 cns->o_viscousStresses, 
                                 cns->o_rhsq);
    } else {
      cns->surfaceKernel(mesh->Nelements, 
                         cns->advSwitch, 
                         mesh->o_sgeo, 
                         mesh->o_LIFTT, 
                         mesh->o_vmapM, 
                         mesh->o_vmapP, 
                         mesh->o_EToB,
                         currentTime, 
                         mesh->o_x, 
                         mesh->o_y,
                         mesh->o_z, 
                         cns->mu,
			 intfx, intfy, intfz,
                         cns->o_rkq, 
                         cns->o_viscousStresses, 
                         cns->o_rhsq);
    }
    
    // update solution using Runge-Kutta
    // rkrhsq_rk = rhsq
    // if rk==6 
    //   q = q + dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
    //   rkerr = dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
    cns->rkUpdateKernel(mesh->Nelements, 
                        rk,
                        mesh->dt, 
                        cns->o_rkA, 
                        cns->o_rkE, 
                        cns->o_q,
                        cns->o_rhsq, 
                        cns->o_rkrhsq, 
                        cns->o_rkq,
                        cns->o_rkerr);
  }
}

void cnsDopriOutputStep(cns_t *cns, const dfloat time, const dfloat dt, const dfloat outTime, occa::memory o_outq){

  mesh_t *mesh = cns->mesh;

  dfloat theta = (outTime-time)/dt; //should have 0<theta<=1

  dfloat *rkB = cns->rkA + 6*cns->Nrk; //the b array is just the last row of A for DOPRI5

  cns->rkoutB[0] = theta*theta*theta*(3-2*theta)*rkB[0] - theta*theta*(theta-1)*(theta-1)*5*((2558722523-31403016*theta)/11282082432) + theta*(theta-1)*(theta-1);
  cns->rkoutB[1] = 0.;
  cns->rkoutB[2] = theta*theta*theta*(3-2*theta)*rkB[2] + theta*theta*(theta-1)*(theta-1)*100*((882725551-15701508*theta)/32700410799);
  cns->rkoutB[3] = theta*theta*theta*(3-2*theta)*rkB[3] - theta*theta*(theta-1)*(theta-1)*25*((2558722523-31403016*theta)/11282082432);
  cns->rkoutB[4] = theta*theta*theta*(3-2*theta)*rkB[4] + theta*theta*(theta-1)*(theta-1)*32805*((23143187-3489224*theta)/199316789632);
  cns->rkoutB[5] = theta*theta*theta*(3-2*theta)*rkB[5] - theta*theta*(theta-1)*(theta-1)*55*((29972135-7076736*theta)/822651844);
  
  cns->rkoutB[6] = theta*theta*(theta-1) + theta*theta*(theta-1)*(theta-1)*10*((7414447-829305*theta)/29380423);


  cns->o_rkoutB.copyFrom(cns->rkoutB);

  cns->rkOutputKernel(mesh->Nelements, 
                      cns->Nrk,
                      mesh->dt, 
                      cns->o_rkoutB, 
                      cns->o_q,
                      cns->o_rkrhsq, 
                      o_outq);
}


void cnsLserkStep(cns_t *cns, setupAide &newOptions, const dfloat time){

  mesh_t *mesh = cns->mesh;
  
  // Low storage explicit Runge Kutta (5 stages, 4th order)
  int advSwitch = 1;//(tstep>100);
    
  for(int rk=0;rk<mesh->Nrk;++rk){
      
    dfloat currentTime = time + mesh->rkc[rk]*mesh->dt;

    dfloat fx, fy, fz, intfx, intfy, intfz;
    cnsBodyForce(currentTime , &fx, &fy, &fz, &intfx, &intfy, &intfz);
    
    // extract q halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nfields;
        
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_q, cns->o_haloBuffer);
        
      // copy extracted halo to HOST 
      cns->o_haloBuffer.copyTo(cns->sendBuffer);      
        
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
    }
      
    // now compute viscous stresses
    cns->stressesVolumeKernel(mesh->Nelements, 
                              mesh->o_vgeo, 
                              mesh->o_Dmatrices, 
                              cns->mu,			      
                              cns->o_q, 
                              cns->o_viscousStresses);
      
    // wait for q halo data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
        
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_q.copyFrom(cns->recvBuffer, cns->haloBytes, offset);
    }
      
    cns->stressesSurfaceKernel(mesh->Nelements, 
                               mesh->o_sgeo, 
                               mesh->o_LIFTT,
                               mesh->o_vmapM, 
                               mesh->o_vmapP, 
                               mesh->o_EToB, 
                               currentTime,
                               mesh->o_x, 
                               mesh->o_y,
                               mesh->o_z, 
                               cns->mu,
			       intfx, intfy, intfz,
                               cns->o_q, 
                               cns->o_viscousStresses);
      
    // extract stresses halo on DEVICE
    if(mesh->totalHaloPairs>0){
      int Nentries = mesh->Np*cns->Nstresses;
          
      mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
        
      // copy extracted halo to HOST 
      cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
          
      // start halo exchange
      meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
    }
      
    // compute volume contribution to DG cns RHS
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {

      cns->cubatureVolumeKernel(mesh->Nelements, 
                                advSwitch,
				fx, fy, fz,
                                mesh->o_vgeo,
                                mesh->o_cubvgeo, 
                                mesh->o_cubDWmatrices,
                                mesh->o_cubInterpT,
                                mesh->o_cubProjectT,
                                cns->o_viscousStresses, 
                                cns->o_q, 
                                cns->o_rhsq);
    } else {
      cns->volumeKernel(mesh->Nelements, 
                        advSwitch,
			fx, fy, fz,
                        mesh->o_vgeo, 
                        mesh->o_Dmatrices,
                        cns->o_viscousStresses, 
                        cns->o_q, 
                        cns->o_rhsq);
    }

    // wait for halo stresses data to arrive
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
        
      // copy halo data to DEVICE
      size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
      cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
    }
      
    // compute surface contribution to DG cns RHS (LIFTT ?)
    if (newOptions.compareArgs("ADVECTION TYPE","CUBATURE")) {
      cns->cubatureSurfaceKernel(mesh->Nelements, 
                                 advSwitch,
                                 mesh->o_vgeo, 
                                 mesh->o_cubsgeo, 
                                 mesh->o_vmapM, 
                                 mesh->o_vmapP, 
                                 mesh->o_EToB,
				 mesh->o_intInterpT,
                                 mesh->o_intLIFTT, 
                                 currentTime, 
                                 mesh->o_intx, 
                                 mesh->o_inty,
                                 mesh->o_intz, 
                                 cns->mu,
				 intfx, intfy, intfz,
                                 cns->o_q, 
                                 cns->o_viscousStresses, 
                                 cns->o_rhsq);
    } else {
      cns->surfaceKernel(mesh->Nelements, 
                         advSwitch, 
                         mesh->o_sgeo, 
                         mesh->o_LIFTT, 
                         mesh->o_vmapM, 
                         mesh->o_vmapP, 
                         mesh->o_EToB,
                         currentTime, 
                         mesh->o_x, 
                         mesh->o_y,
                         mesh->o_z, 
                         cns->mu,
			 intfx, intfy, intfz,
                         cns->o_q, 
                         cns->o_viscousStresses, 
                         cns->o_rhsq);
    }
        
    // update solution using Runge-Kutta
    cns->updateKernel(mesh->Nelements, 
                      mesh->dt, 
                      mesh->rka[rk], 
                      mesh->rkb[rk], 
                      cns->o_rhsq, 
                      cns->o_resq, 
                      cns->o_q);
  }
}
